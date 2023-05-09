import numpy as np
import pandas as pd
import pickle
import copy

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.path import Path
from matplotlib.colors import ListedColormap
import matplotlib.ticker as mtick
import seaborn as sns

import cobra
from cobra import Metabolite, Reaction, Model
from cobra.util.array import create_stoichiometric_matrix
from cobra.util.solver import linear_reaction_coefficients
from cobra.flux_analysis import flux_variability_analysis
from cobra.io.mat import _cell

#from benpy import vlpProblem
#from benpy import solve as bensolve

from scipy.sparse import lil_matrix
from scipy.spatial import Delaunay
from scipy.spatial import distance
import scipy.io as sio

from functools import reduce
from collections import OrderedDict

from eco_utils import *

from sklearn.cluster import DBSCAN, OPTICS,SpectralClustering, AffinityPropagation

from scipy.cluster.hierarchy import fcluster
from scipy.cluster import hierarchy

class Ecosystem:

    def __init__(self, models, prefixes, community_name = "My community model", community_id= "community",
                 pool_bounds = (-1000,1000), keep_members=True, print_report = True, solver='gurobi'):

        self.size = len(models)
        self.name = community_name
        self.prefixes = prefixes
        self.conflicts = None
        self.cmodel = None
        self.exchange_metabolite_info = None
        self.rxn2cluster = None
        self.qual_vector_df = None
        self.fva_results = None
        self.points = None
        self.pfractions = None
        self.feasible_points = None
        
        self.coupled_rxns    = {prefix:None for prefix in self.prefixes}      
        self.member_blocked  = {prefix:None for prefix in self.prefixes}
        self.member_rxns = None
        
        
        
        print("0. Copying member models ...")
        self.models = [m.copy() for m in models]

        
        # --- step 1 ---
        # Store information about exchanged metabolites in each member model
        # information is stored in 'exchange_metabolite_info', 
        # Exchange metabolites with conflicting data between member models (different formulas, charges,
        # or names) are stored in 'conflicts'
        print("1. Storing exchanged metabolites information ...")
        self._get_members_original_exchange_info()
        
        
        # --- step 2 ---
        # Change rxn & met ids and compartments of member models to prefixes[i]:id, prefixes[i]:compartment
        print("2. Updating member model objects ids ...")
        for member_index in range(self.size):
            self._rename_member_model(member_index)
        
        # --- step 3 ---
        # Store member models original objectives {rxn: coef}
        print("3. Storing member objectives ...")  
        self.objectives = list()       
        for ix, m in enumerate(self.models):
            #2.1 objectives
            member_lr_coefs =  linear_reaction_coefficients(m)
            if len(member_lr_coefs) != 1:
                raise RuntimeError("More than one reaction in %s objective function. Not supported!!Set single objective reaction" 
                      % self.prefixes[ix])                
            member_objectives = { r.id:coeff for r,coeff in member_lr_coefs.items()} # should be single key,value dict   
            self.objectives.append(member_objectives)
         
        
        # --- step 4 ---
        # Create new cobrapy Model for the community by merging member models. Use one of the members' objective 
        # function as default objective for the community model
        print("4. Merging member models ...")
        cmodel = Model(id_or_model = community_id, name=community_name)
        for m in self.models:
            cmodel.merge(m)
        
        all_objectives = dict()
        for member in self.objectives:
            print(member)
            for k in member.keys():
                rxn_k = cmodel.reactions.get_by_id(k)
                all_objectives[rxn_k] = member[k] 
            
            
        cmodel.solver = solver
        #print(type(self.objectives[0]))
        #member_objective = self.objectives[0] # we take the objective of one member model      
        #cmodel.objective = { cmodel.reactions.get_by_id(rid) : coeff  for rid, coeff in member_objective.items()}
        cmodel.objective = all_objectives
                
            
        self.cmodel = cmodel       
        
    
        # --- step 5 ---
        # Add pool to community model (including compartment, shared pool mets, pool exchanges
        #                              and update of member exchange reactions)
        # By default pool exchanges are set with bounds (-1000,1000) and individual members' exchange
        # reaction bounds are set to the same values. 
        print("5. Creating pool compartment, metabolites and reactions ...")
        self._add_pool(pool_bounds)
        
        
        # --- step 6 ----
        # Add compartment names. They get lost during merge (step 4)
        print("6. Updating compartment names ...")
        ccomps = self.cmodel.compartments
        for m in self.models:
            for comp in m.compartments:
                ccomps[comp] = m.compartments[comp]
        self.cmodel.compartments = ccomps    
        
        
        # --- step 7 ---
        # Non-blocked reaction list:
        self.non_blocked = None
        print("Done. Community model created")
        
        # --- optional steps ----
        
        # Print summary report if print_report
        if print_report:
            self.print_build_report(pool_bounds)
        # Member models are kept unless keep_members is False:
        if not keep_members:
            self.clean_members()    
    
    def _get_members_original_exchange_info(self):
        
        if self.exchange_metabolite_info is not None:
            print('Exchange information already stored.')     
            return
        
        exchange_metabolite_info = dict()

        #store information of all exchanged metabolites:
        for ix, model in enumerate(self.models):
            prefix = self.prefixes[ix]
            for ex in model.boundary: 
                nmetabolites =len(ex.metabolites)
                if nmetabolites > 1:
                    print("Warning: More than one metabolite in exchange reaction %s from %s model" % (ex.id,prefix))
                    print("Not supported! Use standard of single metabolite exchanges!")
                    return None                
                for m in ex.metabolites:
                    
                    if m.id not in exchange_metabolite_info:
                        exchange_metabolite_info[m.id] = dict()
                    #NJ changed %s:%s to %s_%s    
                    exchange_metabolite_info[m.id][prefix] = {
                        'm_id': "%s_%s" % (prefix,m.id),
                        'name'  :m.name,
                        'formula':m.formula,
                        'charge' :m.charge,
                        'ex_id': "%s_%s" % (prefix,ex.id),
                        'bounds' : ex.bounds,                      
                    } 
        #check for inconsistencies in exchanged metabolites among member models and fill missing
        #values in exchange_metabolite_info
        
        conflicts = list()
        for mid in exchange_metabolite_info:
            names    = set()
            formulas = set()
            charges  = set()
            info = exchange_metabolite_info[mid]
            for prefix in self.prefixes:
                if prefix in info:
                    formulas.add(info[prefix]['formula'])
                    charges.add(info[prefix]['charge'])
                    names.add(info[prefix]['name'])
                else:
                    exchange_metabolite_info[mid][prefix]={
                        'm_id'  : None,
                        'name'  : None,
                        'formula': None,
                        'charge' : None,
                        'ex_id' : None,
                        'bounds' : None,                                       
                    } 
            
            if len(formulas) > 1 or len(charges) > 1 or len(names) > 1:
                conflicts.append(mid) 
                print("Conflicts for exchange metabolite %s:" % mid)
                print("Formulas:")
                print(formulas)
                print("Charges:")
                print(charges)
                print("Names:")
                print(names)
                print("----------")       
            
        self.conflicts = conflicts
        self.exchange_metabolite_info = exchange_metabolite_info
        
        
        return
    
    def print_build_report(self, pool_bounds):
        """ Prints community construction report """
        models = self.models
        cmodel = self.cmodel
        prefixes = self.prefixes
        print("Created community model from %d member models:" % self.size)               
 
        print("General stats:")
        if models is not None:
            for i in range(self.size):
                print("model (%d):" % i) 
                print( "\t id = %s, name = %s , prefix= %s" % (models[i].id,models[i].name, prefixes[i]))
                rxns = len(models[i].reactions)
                comps = len(models[i].compartments)
                mex = len(models[i].boundary)
                print( "\t\t reactions = %d" % rxns)
                print( "\t\t exchange metabolites = %d" % mex)
                print( "\t\t compartments = %d" % comps)
        print("community model:")
        print( "\t id = %s, name = %s " % (cmodel.id, cmodel.name))
        print( "\t\t reactions = %d" % len(cmodel.reactions))
        print( "\t\t exchange metabolites = %d" % len(self.exchange_metabolite_info)) 
        print( "\t\t compartments = %d" % len(cmodel.compartments))        
        print("Exchange metabolite conflicts (formula, charge, name, id) = %d" % len(self.conflicts))   
        print("Community exchange reaction bounds: %s" % str(pool_bounds))
    
    def clean_members(self):
        """ Deletes member individual models"""
        members = self.models
        self.models = None
        del members
        
    def _rename_member_model(self, member_index):
        """ Updates a member model by renaming its metabolite and reaction ids and also 
            its compartment names and ids. 
            The corresponding community member prefixes are used:
                original_rxn_id --> member_prefix:original_rxn_id
                original_met_id --> member_prefix:original_met_id
                original_compartment_id --> member_prefix:original_compartment_id
                original_compartment_name --> member_prefix original_compartment_name

            member_index: index of member model in self.models list. 
         """

        model = self.models[member_index]
        model_prefix = self.prefixes[member_index]
        comp_dict = dict()        

        # Updating metabolite ids and compartment ids and names        
        for m in model.metabolites:
            mid = "%s_%s" % (model_prefix,m.id)#NJ changed %s:%s to %s_%s
            m.id = mid
            mcomp = "%s_%s" % (model_prefix,m.compartment)#NJ changed %s:%s to %s_%s
            if mcomp not in comp_dict:
                
                cname = model.compartments[m.compartment] 
                if len(cname)== 0: # in case compartment doesn't have a name
                    cname = m.compartment
                comp_dict[mcomp]= "%s %s"   % (model_prefix, cname)
            m.compartment = mcomp  
        model.compartments = comp_dict

        #Updating reaction ids
        for rxn in model.reactions:
            rid = "%s_%s" % (model_prefix,rxn.id) #NJ changed %s:%s to %s_%s
            rxn.id = rid
        model.repair()

        return 

    def _add_pool(self,pool_bounds = (-1000,1000)):
        """ 
        Adds pool to community model (including compartment, shared pool mets, pool exchanges
                                      and update of member exchange reactions)
         By default pool exchanges are set with bounds (-1000,1000) 
         Individual member exchange reactions bounds are not modified!      
        """
        cmodel = self.cmodel
        exchange_metabolite_info = self.exchange_metabolite_info
        
        compartments = cmodel.compartments
        compartments['pool']="Community pool"

        #Creating and adding pool metabolites to community model
        pool_mets = list()
        #NJ changed %s:%s to %s_%s in pool
        for mid in exchange_metabolite_info:
            pid = "%s_pool" % mid    
            for prefix in self.prefixes:
                info = exchange_metabolite_info[mid][prefix]
                if info['m_id'] is not None: #info is taking from one of the member models
                        name     = info['name']
                        charge   = info['charge']
                        formula  = info['formula']
                        #print('Adding metabolite %s to pool as %s' % (mid,pid))
                        if mid in self.conflicts:
                            print("Warning: Adding pool metabolite %s with conflicting information in member models" % mid)
                            print("Using data from %s model" % prefix)  
                        mpool = Metabolite(pid, formula = formula, charge = charge, 
                                           name = name, compartment = 'pool')
                        pool_mets.append(mpool)                      
                        break   
        
        cmodel.add_metabolites(pool_mets)             
       
        # Creating and adding pool exchange reactions to community model        

        pool_rxns = list() 
        for m in pool_mets:
            #NJ changed : to _
            mid = m.id.replace('_pool', '')
            ex_rxn = Reaction(
                id   = "EX_%s" % mid,
                name = "Pool %s exchange" % m.name,
                subsystem = "Exchange",
                lower_bound = pool_bounds[0],
                upper_bound = pool_bounds[1]
            )
            ex_rxn.add_metabolites({m: -1.0})        
            pool_rxns.append(ex_rxn)
            
        cmodel.add_reactions(pool_rxns)
        
        # Updating compartment names:
        cmodel.compartments = compartments        
        
        cmodel.repair()
         

        # Updating original exchange reactions from member models:
        # from met_e:model1 <--> to met_e:model1 <--> met_e:pool               
        # Setting bounds to pool bounds
        #NJ changed %s:%s to %s_%s
        for mid in exchange_metabolite_info:
            pid = "%s_pool" % mid 
            mpool = cmodel.metabolites.get_by_id(pid)
            for prefix in exchange_metabolite_info[mid]:
                info = exchange_metabolite_info[mid][prefix]
                if info['m_id'] is not None:
                    mex  = cmodel.metabolites.get_by_id(info['m_id'])
                    ex   = cmodel.reactions.get_by_id(info['ex_id'])
                    coeff = ex.metabolites[mex]
                    ex.add_metabolites({mpool: -coeff})
                
        return         

    def set_pool_bounds(self, pool_metabolites, bioCons = None):
        """ 
        Changes pool exchange reaction bounds of metabolites in pool_metabolites.
        If a metabolite is not part of the pool a warning is raised and nothing is done.
        
        pool_metabolites: dictionary. 
                          keys: metabolite ids (without prefix, i.g., glyc_e) 
                          values: reaction bounds for the corresponding exchange reactions.
        factor: int. Factor to be considered in exchange reactions of the organisms of the community.
                    An organism can consume as maximum as value(from dict)*factor*fraction of that organism
                    in the community. Hence, bigger values make this constraint more flexible
        """
        
        for mid in pool_metabolites:
            if mid in self.exchange_metabolite_info:
                #Changing pool exchange reaction bounds
                rid = "EX_%s" % mid
                ex = self.cmodel.reactions.get_by_id(rid)
                ex.bounds= pool_metabolites[mid]
                #Changing members exchange reaction bounds
                for member in self.exchange_metabolite_info[mid]:
                    ex_id = self.exchange_metabolite_info[mid][member]['ex_id']

                    if ex_id is not None:
                        ex = self.cmodel.reactions.get_by_id(ex_id)
                        nlb = pool_metabolites[mid][0]
                        if nlb <= 0 and bioCons is not None:
                            #ex.lower_bound = nlb
                            ex.lower_bound = bioCons
                
                
                
            else:
                print("Skipping %s ..." % mid)
                print("Warning: %s is not a part of pool metabolites." % mid)
            self.non_blocked = None    

    def set_member_exchange_bounds(self, member_prefix, exchange_metabolites):
        """ 
        Changes member exchange reaction bounds of metabolites in exchange_metabolites.
        If a metabolite is not part of the exchanges a warning is raised and nothing is done.
        
        member_prefix: member model whose exchange reactions are modified
        exchange_metabolites: dictionary. 
                          keys: metabolite ids (without prefix, i.g., glyc_e)
                          values: reaction bounds for the corresponding exchange reactions.
        """            
            
            
        df =  self.get_exchange_df('ex_id') #id of member exchange reactions   
        
        for m in exchange_metabolites:
            new_bounds = exchange_metabolites[m]
            if m in df.index:
                rid = df.loc[m, member_prefix]
                if rid is not None:
                    rxn = self.cmodel.reactions.get_by_id(rid)
                    rxn.bounds = new_bounds
                    self.exchange_metabolite_info[m][member_prefix]['bounds'] = new_bounds
                else:
                    print("No exchange reaction for %s in %s. Skypping..." % (m, member_prefix))
            else:
                print("No exchange or pool reactions for %s. Skypping" % m)
                 
    def show_member_exchanges(self, mids=None):
        df = self.get_exchange_df('bounds')
        if mids is not None:
            df = df.loc[mids]
        return df    
    
    def get_exchange_df(self, k): 
    
        index = sorted(self.exchange_metabolite_info.keys())
        columns = sorted(self.prefixes)
        rows=list()
        for m in index:
            info = self.exchange_metabolite_info[m]
            row = [info[prefix][k] for prefix in columns]
            rows.append(row)
        
        df = pd.DataFrame(data=rows, index=index, columns=columns)
        return df
    
    def _to_vlp(self,**kwargs):        
        """Returns a vlp problem from EcosystemModel"""
        # We are using bensolve-2.0.1:
        # B is coefficient matrix
        # P is objective Matrix
        # a is lower bounds for B
        # b is upper bounds for B
        # l is lower bounds of variables
        # s is upper bounds of variables
        # opt_dir is direction: 1 min, -1 max
        # Y,Z and c are part of cone definition. If empty => MOLP
        
        cmodel = self.cmodel
        Ssigma = create_stoichiometric_matrix(cmodel, array_type="lil")
        
        vlp = vlpProblem(**kwargs)
        m, n = Ssigma.shape # mets, reactions
        q = self.size # number of members 
        vlp.B = Ssigma
        vlp.a = np.zeros((1, m))[0]
        vlp.b = np.zeros((1, m))[0]
        vlp.l = [r.lower_bound for r in cmodel.reactions] 
        vlp.s = [r.upper_bound for r in cmodel.reactions] 
        
        
        vlp.P = lil_matrix((q, n))
        vlp.opt_dir = -1
        
        for i, member_objectives in enumerate(self.objectives):
            for rid, coeff in member_objectives.items():
                rindex = cmodel.reactions.index(rid)
                vlp.P[i,rindex] = coeff 
                
        vlp.Y = None
        vlp.Z = None
        vlp.c = None
        return vlp  
    
    def mo_fba(self, bensolve_opts = None):
       
        if bensolve_opts is None:
            bensolve_opts = vlpProblem().default_options
            bensolve_opts['message_level'] = 0
        
        vlp_eco = self._to_vlp(options = bensolve_opts)    
        self.mo_fba_sol = bensolve(vlp_eco)
         
    def get_polytope_vertex(self, expand=True):
  
        """
        polytope: pareto front + axes segments + extra segments perpendicular to axes dimensions where 
        pareto solutions don't reach 0 values. 
        (assumption: objective functions can only take positive values)
        """
        
        #1. Front vertex:
        vv = self.mo_fba_sol.Primal.vertex_value[np.array(self.mo_fba_sol.Primal.vertex_type)==1]
        
        n_neg_vals = np.sum(vv<0)
        if n_neg_vals > 0:
            print('warning: Negative values in Pareto Front..')
            print(vv[vv<0])
            print("Changing negative values to zero...")
            vv[vv<0] = 0        
        
        
        #2. origin
        ov = np.zeros((1,self.size))
        
        
        if expand == True:
            #3. Extra points that close polytope (only if points (0,0,0,...,xi_max,0,...0) are not pareto front 
            # points but they are feasible) 
            
            #MP: si hay que incluir estos puntos significa que hay miembros que son givers: i.e. pueden crecer
            # a su máxima tasa y aguantar que otros miembros crezcan también
            # si un punto (0,0,0,...,xi_max,0,...0) no es factible entonces el miembro i necesita que otros crezcan 
            # para poder crecer (needy).
            
            #3.1 Check if points  (0,0,0,...,xi_max,0,...0) are part of the pareto front
            n = self.size - 1
            all_zeros_but_one = np.argwhere(np.sum(vv == 0,axis=1)==n) # index of (0,0,...,xi_max,0,0...0) points
            all_zeros_but_one = all_zeros_but_one.flatten()
    
            # indexes i of non-zero member in (0,0,...,xi_max,0,0...0) pareto points, 
            # i.e. members that are not givers nor needy. 
            non_zero_dims =  np.argmax(vv[all_zeros_but_one,:], axis = 1) 
            
            # givers and/or needy members:
            givers_or_needy_indexes = np.setdiff1d(np.array(range(self.size)), non_zero_dims) 
            gn_total= len(givers_or_needy_indexes)    
        
            #3.2 Check if non-pareto points (0,0,0,...,xi_max,0,...0) are feasible 
            if gn_total >0:
                # max values for giver_or_needy members:
                max_vals = np.max(vv, axis=0)
                cpoints = np.diag(max_vals)
                to_check = cpoints[givers_or_needy_indexes,:]
                
                are_feasible = self.check_feasible(to_check)
                
                ev = to_check[are_feasible,:] 
                polytope_vertex = np.concatenate((vv,ov,ev), axis=0)
            else: 
                polytope_vertex = np.concatenate((vv,ov), axis=0)
        
        else:
                polytope_vertex = np.concatenate((vv,ov), axis=0)
        return polytope_vertex 
    
    def check_feasible(self, point_array, pfraction_array=None, update_bounds=False, **kwargs):
        if not update_bounds:
            r = [self._analyze_point((p,None), analysis = 'feasibility', update_bounds=False, **kwargs) for p in point_array] 
            
        else:
            npoints = point_array.shape[0]
            nfrac   = pfraction_array.shape[0]
            if pfraction_array is None or npoints != nfrac:
                raise RuntimeError("Missing or incomplete member fractions array. Cannot update rxn bounds!!") 
            else:
                coord_frac = [(point_array[ix],pfraction_array[ix]) for ix in range(npoints)]    
                r = [self._analyze_point(x, analysis = 'feasibility', update_bounds = True, **kwargs) for x in coord_frac] 
            
        return r
    
    def calculate_qual_vectors(self,point_array, pfraction_array=None, update_bounds=False, **kwargs):
        
        # Check for reactions selected for FVA and clustering
        if self.rxn2cluster is None:
            print("No reactions previously selected for FVA and clustering!")
            print("Setting reactions to cluster...")
            self.set_cluster_reactions()        
        
        
        
        if not update_bounds:
            print("Warning:Calculating qualitative vectors without updating reaction bounds!")
            r = [self._analyze_point((p,None), analysis = 'qual_fva', update_bounds=False, **kwargs) for p in point_array] 
        else:
            npoints = point_array.shape[0]
            nfrac   = pfraction_array.shape[0]
            
            if pfraction_array is None or npoints != nfrac:
                raise RuntimeError("Missing or incomplete member fractions array. Cannot update rxn bounds!!") 
            else:
                coord_frac = [(point_array[ix],pfraction_array[ix]) for ix in range(npoints)]  
                #list of tuples
                r = [self._analyze_point(x, analysis = 'qual_fva', update_bounds = True, **kwargs) for x in coord_frac] 
         
        return r
    
    def _analyze_point(self, ptuple, analysis = 'feasibility', update_bounds = False, delta = 1e-9):
        
        # ptuple: tuple of two arrays:   (0) Grid point coordinates; 
        #                                (1) grid point member fractions; grid point
        # analysis: type of analysis to run. Options:
        #           'feasibility': check if grid point is feasible
        #           'qual_fva'   : get grid point vector of rxn qualitative values. 
        #                           
                    
        # update_bounds: if True update reaction bounds considering member community fractions 
        #                before analysis 
        # delta: threshold value to consider flux differences as zero when comparing fva min and max values
        #        ('qual_fva' option only)  
        # Returns  
        #  boolean indicating point feasibility ('feasibility' analysis) or  ('qual_val' analysis)
        #  a tuple where the first element is a list of rxns qualitative values for the analyzed grid point
        #  and the second is an array with the corresponding fva results
        def qual_translate(x, delta = 1e-4):
            fvaMax = x.maximum
            fvaMin = x.minimum
           
            
            #3 fixed on positive value
             #-3 fixed on negative value
             
            # -2: max and min <0
            if fvaMax < -delta and fvaMin < -delta:
                ret = -2
                if abs(fvaMax-fvaMin) < delta:
                    ret = -3
            # 0: max and min == 0
            elif (fvaMax >= -delta and fvaMax <=delta) and (fvaMin >=-delta and fvaMin <=delta):
                ret = 0
            # -1: max and min <=0
            elif (fvaMax>= -delta and fvaMax < delta) and (fvaMin < -delta):
                ret = -1
            # 1: max and min >= 0
            elif (fvaMin >=-delta and fvaMin <=delta) and fvaMax >delta:
                ret = 1
            # 2: max and min >0
            elif (fvaMax >delta and fvaMin > delta):
                ret = 2
                if abs(fvaMax - fvaMin) < delta:
                    ret = 3
            elif (fvaMin < -delta and fvaMax > delta):
                ret = 4 #reversible
            else:
                ret = 5 #nans            

            
            ## category - (-3): max = min < 0
            #if (fvaMax - fvaMin) <= delta and fvaMax < -delta:
            #    ret = -3
            
            ## category + (5): max = min > 0
            #elif (fvaMax - fvaMin) <= delta and fvaMin > delta:
            #    ret = 5     
            
            ## category -- (-2): max and min <0     
            #elif fvaMax < -delta and fvaMin < -delta:
            #    ret = -2
            
            ## category 0  (0): max and min == 0
            #elif (fvaMax >= -delta and fvaMax <=delta) and (fvaMin >=-delta and fvaMin <=delta):
            #    ret = 0
            
            ## category -0 (-1): max and min <=0
            #elif (fvaMax>= -delta and fvaMax < delta) and (fvaMin < -delta):
            #    ret = -1
            
            ## category 0+ (1): max and min >= 0
            #elif (fvaMin >=-delta and fvaMin <=delta) and fvaMax >delta:
            #    ret = 1
            
            ## category ++ (2): max and min >0
            #elif (fvaMax >delta and fvaMin > delta):
            #    ret = 2
            
            ## category -+ (3): max > 0 and min < 0     
            #elif (fvaMin < -delta and fvaMax > delta): 
            #    ret = 3            
            
            ## category * (4): something weird happened!
            #else:
            #    ret = 4 #nans            
            return ret          
        
        # 0,0+,0-,++,--,+,-,+-,*
        
        cmodel = self.cmodel
        point = ptuple[0] #grid point coordinates
        print(point)
        point = [point[0]*point[1], (1-point[0])*point[1]]#equivalent of old point 
        pfrac = ptuple[1] #grid point member fraction
            
        out = None
        
        with cmodel:
           # update member reactions bounds if required:
           if update_bounds: 
                print('updating reaction bounds ...')
                # In the case that all objectives are zero, all member fractions are nan
                # Here we set fractions to zero to avoid errors setting bounds
                if np.all(np.isnan(pfrac)): 
                    pfrac = np.array([0.0]*pfrac.size)                
                    
                # reactions are assign to each community member
                self.get_member_reactions()
                for i, member in enumerate(self.prefixes):
                    mfrac = pfrac[i]
                    mrxns = self.member_rxns[member]
           
                    # rxn bounds are updated, accounting for members fractions in the community 
                    
                    for rid in mrxns:
                        r = cmodel.reactions.get_by_id(rid)
                        old_bounds = r.bounds
                        r.bounds =(old_bounds[0] * mfrac, old_bounds[1] * mfrac)

           # fix member objectives to grid point value:
           for ix, member_objectives in enumerate(self.objectives):    
                if len(member_objectives) != 1:
                    raise RuntimeError("Warning: More than one reaction in %s objective function. Not supported!!" 
                      % self.prefixes[ix])        
                #new bounds for member ix objective function reaction:
                new_bounds = (point[ix],point[ix])    
                #newGrid NJ
                #print(point[ix])

            
                #change bounds for each objective reaction
                for rid in member_objectives.keys(): # member_objectives should be single key dictionary
                    rxn = cmodel.reactions.get_by_id(rid)
                    rxn.bounds = new_bounds

                    #set one of objective reactions as community objective
                    #cmodel.objective = rid #commented since the model comes with an objective
            
           
           if analysis == 'feasibility': 
                #check point feasibility
                error_value = -1000
                ob = cmodel.slim_optimize(error_value = error_value)  
                if ob == error_value:
                    out = False
                    print('unfeasible point')
                else:
                    out = True
        
           elif analysis == 'qual_fva':  # here we assume the point is feasible      

                if self.rxn2cluster is None:
                    raise RuntimeError('No reactions selected for fva and clustering!')
                    
                print("running FVA on grid point...")
                print(ptuple)
                
                rxn_fva = flux_variability_analysis(cmodel,reaction_list= self.rxn2cluster)
                
                rxn_fva = rxn_fva.loc[self.rxn2cluster,:] # just to make sure reactions are in the 
                                                         # same order as rxn2cluster
                    
                print("translating to qualitative vector..")
                out = (list(rxn_fva.apply(qual_translate, axis = 1, delta = delta)), #qual_vector 
                             rxn_fva.values) # array with FVA results  
        return out      
        
    def analyze_grid(self, analysis = 'feasibility', update_bounds=True, **kwargs):
        #analysis: type of analysis to run on full grid:
        #   Options:         
        #       feasibility = checks if each grid point is feasible considering member fractions
        #       qual_fva = calculates qualitative vectors for each point in the grid. If feasible 
        #                  points are stored, analysis is run on those points only.   
        #                  FVA results are also stored.   
        
        #step 1: # calculate community distribution (member fractions) for each grid point 
                 # if they are not stored
        if self.pfractions is None:
            self.get_points_distribution()
        
        point_array     = self.points               
        pfraction_array =self.pfractions
        
        #Option 'feasibility':
        if analysis == 'feasibility':
                
            #step 2: run feasibility analysis for all grid points
            npoints = point_array.shape[0]       
            self.feasible_points = self.check_feasible(point_array, pfraction_array, 
                                                       update_bounds=update_bounds, **kwargs)
            
            npoints = point_array.shape[0]
            nfeasible = sum(self.feasible_points)
            print("grid feasible points: %d / %d" % (nfeasible, npoints))
            
        elif analysis == 'qual_fva':
                
            #step 2: check for previously calculated feasible points  
            feasible_points = self.feasible_points
            if feasible_points is None:
                print("Warning: Feasible points have not been calculated. Running qualitative fva over full grid")
                df_index = np.arange(point_array.shape[0])
            else:
                print("Running qualitative fva over grid feasible points...")
                point_array = point_array[feasible_points,:]    
                pfraction_array = pfraction_array[feasible_points,:]     
                df_index =  np.where(feasible_points)[0]
        
            rtuples = self.calculate_qual_vectors(point_array,pfraction_array,update_bounds=update_bounds, **kwargs)
            
            qual_vector_list, fva_results =  map(list, zip(*rtuples))    
            self.qual_vector_df = pd.DataFrame(np.array(qual_vector_list),columns = self.rxn2cluster, index=df_index)
            
            fva_results = np.dstack(fva_results)
            fva_results = np.rollaxis(fva_results,-1)
            
            self.fva_results = fva_results    
                
    def change_reaction_bounds(self, rid, new_bounds):
        cmodel = self.cmodel
        rxn = cmodel.reactions.get_by_id(rid)
        old_bounds = rxn.bounds
        rxn.bounds = new_bounds
        return old_bounds
                
    def build_grid(self,expand = True, numPoints=10, drop_zero = True, ignore_maint = True):
                
        #compute max_com by relaxing constraints such as ATPM
        with self.cmodel:
            if ignore_maint:
                for rxn in self.cmodel.reactions:
                    if rxn.lower_bound > 0:
                        rxn.lower_bound = 0
            
            max_com = self.cmodel.slim_optimize()
        maxs = [1 ,max_com]
        print('Maximum community:'+str(max_com))
        mins = [0,0] #hardcoded 2D 
        size = self.size
        
        #Modify this to have a matrix of nxn points rather than a step (using com_growth and fraction as axis)
        #alternative: define a different step based on points to have
        slices = [np.linspace(mins[i], maxs[i], numPoints) for i in range(size)]
        #slices = [slice(mins[i],maxs[i],step) for i in range(size)]
        #rgrid = np.mgrid[slices]

        print(slices[0])
        print(slices[1])
        rgrid = np.array(np.meshgrid(slices[0], slices[1]))#.T.reshape() #NJ
        print(rgrid)
        #rgrid2columns = [rgrid[i,:].ravel() for i in range(size)]
        rgrid2columns = [rgrid[i,:].ravel() for i in range(size)]
        # array of grid points (x,y,z,...)
        positions = np.column_stack(rgrid2columns)
        #polytope intersection    
        #inside = hull.find_simplex(positions)>=0
        
        # storing grid points inside polytope
        #points = positions[inside] 
        points = positions
        limits = (mins,maxs)
        if drop_zero:
            points = points[1:]
            mins   = np.min(points, axis=0)
            maxs   = np.max(points, axis=0)
        
        
        
        self.points  = points
        #self.step    = step
        self.limits  = (mins,maxs)

    def get_selected_points_distribution(self, prange = None):
        
        if self.points is None:
            raise RuntimeError('Grid points are not set yet!')
            
        
        points = self.points
        #NJ new grid
        #points[0]: f_i
        #points[1]: com_u
        
        if prange is not None:
            points = points[prange]
        
        #com_mu = np.sum(points,axis =1)
        #pfractions = np.apply_along_axis(lambda a, b : a/b, 0, points, com_mu)     
        pfractions = np.array([[p[0], 1-p[0]] for p in points]) #assuming two organisms
        return pfractions        
        
    def get_points_distribution(self):
        pfractions = self.get_selected_points_distribution(prange = None)
        self.pfractions = pfractions

    def calculate_community_growth(self, feasible=False): #NJ DELETE THIS FUNCTION
        cgrowth = np.sum(self.points,axis=1) 
        if feasible:
            if self.feasible_points is None:
                print('feasible points have not been previously established! Returning values for all points')
            else:
                cgrowth =cgrowth[self.feasible_points]
                       
        return cgrowth 
     
    def get_member_reactions(self):
        member_rxns = {x:[] for x in self.prefixes}
        for r in self.cmodel.reactions:
            for member in self.prefixes:
                if r.id.startswith(member):
                    member_rxns[member].append(r.id)
                    break
        self.member_rxns = member_rxns                  
    
    def get_2D_slice(self, prefixes, fixed_values): #NJ DELETE THIS FUNCTION
        if self.size - len(prefixes) != 2:
            raise RuntimeError("Two members with non-fixed values required! No more, no less.")
        
        
        members_to_fix = range(len(prefixes))
        #valores mas cercanos en la grilla a fixed_values
        closest_values = [self._get_closest_grid_value(prefixes[i],fixed_values[i])for i in members_to_fix]
        fixed_member_indexes = [self.prefixes.index(x) for x in prefixes]
        
        
        grid_points = self.points
        #grid_points_values = self.points_values # aqui va matriz de puntos x reacciones con valores cualitativos
                                         # calculados a partir de fva. 
        
        #indices de puntos en grilla donde el valor una dimension es igual a su closest_values
        filtered_indexes_list =  [np.where(grid_points[:,fixed_member_indexes[i]] == closest_values[i])[0] for i in members_to_fix]
        
        #interseccion de indices, i.e., indices slice
        slice_indexes =  reduce(np.intersect1d, filtered_indexes_list)
        #slice 2D de la grilla
        #filtered_points = grid_points[slice_indexes,:] 
        #filtered_values = grid_points_values[slice_indexes,:]
        
        free_members = list(set(self.prefixes).difference(set(prefixes))) # prefixes of members with non-fixed objective values
        free_member_indexes = [self.prefixes.index(x) for x in free_members]
        #Puntos se reducen a las dimensiones de los free members:
        #slice_points = filtered_points[:,free_member_indexes]
        #slice_points_values = filtered_values[:,free_member_indexes]
          
        #return slice_points, slice_points_values , free_member_indexes
        return [slice_indexes, free_member_indexes]     

    def plot_2D_slice(self, prefixes=[], fixed_values=[], parent_cmap='tab20',s=8, xlabel = None, ylabel = None, 
                      figsize=(11,12), to_plot = None, show_edge=False,frac_prefix= None, saveFile='', fractions = None):
        
        """ Plots a 2D slice of the grid. The remaining dimensions of the grid are kept at fixed values.
            Individual organism models of fixed dimensions are those determined by 'prefixes'. Their 
            objective function values are fixed to 'fixed_values'.   
            Points previously found unfeasible are not shown.
        
            prefixes: List of individual organism prefixes to be set at fixed objective function values 
            fixed_values: Fixed values for objective functions of 'prefixes' individual models
            parent_cmap: color map used to show clusters
            to_plot: How to color grid points:
                 'cluster':  Points are colored according to their clusters. Clusters must be previously calculated.  
                             if point feasibility has been run, only feasible points are shown.  
                 'feasible': Points are colored as feasible or unfeasible. Point feasibility must 
                             be previously calculated.
                             
                  'community_growth': Points are colored according to their corresponding community growth rate.
                  'member_fraction' : The community fraction of a particular member (frac_prefix) is used to color grid points
                  
                  None: no particular coloring over grid points. 
            
            xlabel: x label. If None member prefix is used
            ylabel: y label. If None member prefix is used
            
            s: marker size
            figsize: figure size
            show_edge: Draw grid edge.
            
            Returns:       
            slice_points: grid points in slice.
        """    
        
        # get full 2D slice:   
        if self.size == 2:
            if len(prefixes)>0:
                print("Only two members in community!!") 
                print("Full grid will be plotted and fixed values for %s will be ignored..." % str(prefixes))
            free_member_prefixes = self.prefixes
            free_member_indexes = [0,1]
            full_slice_indexes = np.arange(len(self.points)) 
            
        else:              
            full_slice_indexes, free_member_indexes = self.get_2D_slice(prefixes, fixed_values)
            free_member_prefixes = [self.prefixes[x] for x in free_member_indexes]
        
        
        # get edge of full slice
        lines = None
        full_slice_points = self.points[full_slice_indexes,:][:,free_member_indexes]
        if show_edge:
            full_slice_points = self.points[full_slice_indexes,:][:,free_member_indexes]
            lines = self._draw_slice_edge(full_slice_points)        

        # get slice feasible points and their member distribution:
        if self.feasible_points is not None: 
            feasible_indexes    = np.where(self.feasible_points)[0]
            feasible_points     = self.points[feasible_indexes]
            feasible_pfractions = self.pfractions[feasible_indexes]
            
            slice_indexes    = np.isin(feasible_indexes, full_slice_indexes)
            slice_points     = feasible_points[slice_indexes,:][:,free_member_indexes]  
            #slice_pfractions = feasible_pfractions[slice_indexes]
            slice_pfractions = feasible_pfractions
        else:
            slice_indexes    = full_slice_indexes
            slice_points     = self.points[full_slice_indexes,:][:,free_member_indexes] 
            slice_pfractions = self.pfractions[full_slice_indexes]

        #axes labels    
        if xlabel is None:
            xlabel = free_member_prefixes[0]
        if ylabel is None:    
            ylabel = free_member_prefixes[1] 

        #Different plot coloring cases  
        
        # 1. Nothing is colored
        if to_plot is None: 
        
             self._no_coloring_plot(slice_points, lines=lines, xlabel = xlabel, ylabel= ylabel, figsize=figsize, s=s, saveFile=saveFile)
             
        # 2. Points are colored as feasible and unfeasible        
        elif to_plot == 'feasible':
            if  self.feasible_points is None:
                raise RuntimeError('Points feasibility analysis has not been run! Required for plot!')
                 
            aux = np.array(self.feasible_points)[full_slice_indexes]
            slice_colors = aux + 1
            k = 2
            color_labels = ['','unfeasible', 'feasible']
            
            self._categorical_coloring_plot(full_slice_points, slice_colors, parent_cmap, k, color_labels, lines=lines, xlabel = xlabel, ylabel= ylabel,
                                   figsize=figsize, s=s,shrink=0.5)
            #self._categorical_coloring_plot(slice_points, slice_colors, parent_cmap, k, color_labels, lines=lines, xlabel = xlabel, ylabel= ylabel,
            #                       figsize=figsize, s=s,shrink=0.5)
            
        # 3. Points are colored according to their community growth values    
        elif to_plot == 'community_growth':
        
            cgrowth = self.calculate_community_growth(feasible=True)
            slice_colors = cgrowth[slice_indexes]
            
                        
            self._gradient_coloring_plot(slice_points, slice_colors, parent_cmap, color_label = 'community growth', 
                                         xlabel = xlabel, ylabel= ylabel, figsize=figsize, s=s, shrink=0.5, lines = lines)    
            
        # 4. Points are colored according to their clusters 
        elif to_plot == 'cluster':

            slice_colors = self.clusters[slice_indexes] #slice points clusters
            k = self.k            
            color_labels = ['']+['c'+ str(x+1) for x in range(k)] 
            pfrac = self.pfractions
            if fractions is not None:
                fracBool = np.isclose(pfrac, fractions)
                print(pfrac[fracBool])
                print(sum(fracBool))
                #fracBool = pfrac == fractions
                frac = [f[0] and f[1] for f in fracBool]
                pointsF = [[self.points[i][0],self.points[i][1]] for i,f in enumerate(frac) if f]
                print(pointsF)
            else:
                pointsF = None

            self._categorical_coloring_plot(slice_points, slice_colors, parent_cmap, k, color_labels, lines=lines, 
                                            xlabel = xlabel, ylabel= ylabel,figsize=figsize, s=s,shrink=1, saveFile = saveFile, pointsF=pointsF)
            
            
             
        # 5. Points are colored according to a selected member community fraction
        elif to_plot == 'member_fraction':
            if frac_prefix is None or frac_prefix not in self.prefixes:
                 raise RuntimeError('Missing valid member prefix: frac_prefix')
             
            member_index = self.prefixes.index(frac_prefix)
            slice_colors = slice_pfractions[:,member_index]
            color_label = "%s fraction" % frac_prefix
            
            self._gradient_coloring_plot(slice_points, slice_colors, parent_cmap, color_label = color_label, 
                                         xlabel = xlabel, ylabel= ylabel, figsize=figsize, s=s, shrink=0.5, lines = lines)     
    
    @staticmethod
    def _draw_slice_edge(full_slice_points):
        #Delaunay triangulation para slice
        slice_delaunay = Delaunay(full_slice_points)        
        
        #Puntos en el borde de slice:
        edges = set()
        edge_points = []   
        for i, j in slice_delaunay.convex_hull:
            if (i, j) not in edges and (j, i) not in edges:
                edges.add( (i, j) )
                edge_points.append(slice_delaunay.points[ [i, j] ])
        

        #plot de borde de slice
        lines = LineCollection(edge_points, color='dimgrey')    
        return lines  
                 
    @staticmethod
    def _no_coloring_plot(slice_points, lines=None, xlabel = None, ylabel= None, figsize=(11,12), s=16, saveFile=''):
        
        myfig = plt.figure(figsize=figsize)            
        if lines is not None:
            plt.gca().add_collection(lines)    
            
        ps = plt.scatter(slice_points[:,0],slice_points[:,1],s=s) 
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if len(saveFile)>0:
            plt.savefig(str(saveFile)+'.pdf')
        plt.show()        
                 
    @staticmethod
    def _categorical_coloring_plot(slice_points, slice_colors, parent_cmap, k, color_labels,  lines = None, 
                                   xlabel = None, ylabel= None, figsize=(11,12), s=16,shrink=1, saveFile='', pointsF=None):
        myfig = plt.figure(figsize=figsize)            
        if lines is not None: 
            plt.gca().add_collection(lines)
        cmap1=plt.get_cmap(parent_cmap,k)
        #np.random.shuffle(cmap1.colors) #colors are not random
        vmin=0.5
        vmax =k + 0.5 
        #ps = plt.scatter(slice_points[:,0],slice_points[:,1], c = slice_colors,cmap=cmap1, s=s,vmin=vmin,vmax=vmax)

        if pointsF is not None:
            x = [0]+[p[0] for p in pointsF]
            y = [0]+[p[1] for p in pointsF]

            ps = plt.scatter(slice_points[:,0],slice_points[:,1], c = slice_colors,cmap=cmap1, s=s,vmin=vmin,vmax=vmax, alpha=0.5)

            plt.plot(x,y,'k')
        else:
            ps = plt.scatter(slice_points[:,0],slice_points[:,1], c = slice_colors,cmap=cmap1, s=s,vmin=vmin,vmax=vmax)
       
        cbar = myfig.colorbar(ps, ticks=np.arange(k+1),shrink=shrink)        
        cbar.ax.set_yticklabels(color_labels)  
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if len(saveFile)>0:
            plt.savefig(str(saveFile)+'.png')
        plt.show()           

    
    
    @staticmethod
    def _gradient_coloring_plot(slice_points, slice_colors, parent_cmap, color_label = None, xlabel = None, ylabel= None,
                                   figsize=(11,12), s=16, shrink=1, lines = None):
        myfig = plt.figure(figsize=figsize)            
        if lines is not None:
            plt.gca().add_collection(lines)
    
        ps = plt.scatter(slice_points[:,0],slice_points[:,1], c = slice_colors,cmap=parent_cmap,s=s)            
        cbar = myfig.colorbar(ps, shrink=shrink)
        #cbar.ax.set_ylabel(color_label, rotation=270)
        cbar.set_label(color_label, rotation=270, labelpad= 20)#, labelpad=0.5)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.show()         
    
    def _get_closest_grid_value(self,prefix,fixed_value):
        member_index = self.prefixes.index(prefix)
        member_min = self.limits[0][member_index]
        member_max = self.limits[1][member_index]
        
        if fixed_value < member_min or fixed_value > member_max:
            raise RuntimeError("Value %d for %s out of range" % (fixed_value, prefix))
        
        
        shifted_value = fixed_value - member_min
        n_steps = shifted_value//self.step          
        
        p1 = n_steps * self.step
        p2 = (n_steps + 1) * self.step  
        
        if (shifted_value - p1) < (p2 - shifted_value):
            closest_value = member_min + p1
        else:
            closest_value = member_min + p2        
        
        return closest_value    
    
    def get_blocked_reactions(self):
        blocked = cobra.flux_analysis.find_blocked_reactions(self.cmodel)
        return blocked     
        
    def get_non_blocked_reactions(self):
        blocked = cobra.flux_analysis.find_blocked_reactions(self.cmodel)
        all_ids = [x.id for x in self.cmodel.reactions]
        non_blocked = set(all_ids).difference(set(blocked))
        self.non_blocked = non_blocked
        
    def clusterPoints(self, method, numeric_delta=1e-4, vector = 'qual_vector', run_fva = True , **kwargs):
        cmodel = self.cmodel

        if self.qual_vector_df is None:
            print("No qualitative FVA values stored. Run qual_fva analysis first!")
            return

        print("Calculating jaccard distances between grid points...") 
        
        #NJT to use qual_vector as well as bin_vector if required
        if vector =='qual_vector':
            z = self.qual_vector_df.values 
        elif vector == 'bin_vector' and self.bin_vector_df is not None:
            z = self.bin_vector_df.values
        
        distance_metric = 'jaccard'
        dvector = distance.pdist(z,distance_metric)
        dmatrix = distance.squareform(dvector)
        
        print("Clustering grid points ...") 
        # Clustering methods
        # hierarchical + fcluster
        # dbscan (eps,min_samples)
        # optics
        # AffinityPropagation
        # SpectralClustering 
       
        if method == 'hierarchical':
            k, clusters = getHIERARCHICALclusters(dvector,**kwargs)
        elif method == 'dbscan':
            k, clusters = getDBSCANclusters(dmatrix,**kwargs)  
        elif method == 'optics':
            k, clusters = getOPTICSclusters(dmatrix,**kwargs)
        elif method == 'SpectralClustering':
            k, clusters = getSCclusters(dmatrix,**kwargs)          
        elif method == 'AffinityPropagation':
            k, clusters = getAPclusters(dmatrix, **kwargs)
            
        self.k = k
        self.clusters = clusters    
        print("Done!")        
        
    def write_fca_input(self,prefix,file_dir, discard_zero_bound=True):
     
        file_name = "%s/%s_fca_input.mat" % (file_dir,prefix)
        cmodel = self.cmodel 
            
        rxns = cmodel.reactions
        mets = cmodel.metabolites
        stoich_mat = create_stoichiometric_matrix(cmodel)
        rids = np.array(rxns.list_attr("id"))
        mids = np.array(mets.list_attr("id"))
        rev  = np.array(rxns.list_attr("reversibility"))*1
    
        #discard reactions from pool and other members and also reactions from prefix with zero lower and upper bounds 
        #these last reactions are added to blocked.
        
        to_discard = []
        blocked = [] 
                
        for ix in range(len(rids)):
            rid = rids[ix]
            if not rid.startswith(prefix):
                to_discard.append(ix)
            else:        
                if discard_zero_bound:            
                    r = cmodel.reactions.get_by_id(rid)
                    if r.lower_bound == 0 and r.upper_bound == 0 :
                        to_discard.append(ix)
                        blocked.append(rid)     
        
        if len(to_discard)>0:
            rids = np.delete(rids,to_discard)
            rev = np.delete(rev,to_discard)
            stoich_mat = np.delete(stoich_mat, to_discard, axis=1)
                  
        
        #discard metabolites from pool and other members:
        to_discard = []
        for ix in range(len(mids)):
            mid=mids[ix]
            if not mid.startswith(prefix):
                to_discard.append(ix)
        
        if len(to_discard)>0:  
            mids = np.delete(mids,to_discard)
            stoich_mat = np.delete(stoich_mat, to_discard, axis=0)    
        
        rids =list(rids)
        mids =list(mids)
        
        
        #create mat objects for FCA and stored them in an output file 
        mat  = OrderedDict()    
        mat["Metabolites"] = _cell(mids)
        mat["Reactions"] = _cell(rids)
        mat["stoichiometricMatrix"] = stoich_mat
        mat["reversibilityVector"] = rev 
        
        mat2 = OrderedDict()
        #mat2['bound_blocked'] = _cell(blocked)
        varname1 = 'fca_input'
        varname2 = 'bound_blocked'
        #varname2 = "%s_bound_blocked" % prefix
        #sio.savemat(file_name, {varname1: mat, varname2:mat2}, oned_as="column")
        sio.savemat(file_name, {varname1: mat, varname2:_cell(blocked)}, oned_as="column")
        print("%s FCA's input in %s" % (prefix, file_name))
        print("   stoichiometric matrix : %s" % str(stoich_mat.shape))
        
    def store_fca_results(self,prefix,fca_file):

        mat_contents = sio.loadmat(fca_file)
        fctable = mat_contents['fctable']
        blocked = mat_contents['blocked'][0]
        rxn_ids = np.array([rid[0][0] for rid in mat_contents['rxns']])
        bound_blocked = np.array([rid[0][0] for rid in mat_contents['bound_blocked']])
        blocked_ids = rxn_ids[blocked==1]
        blocked_ids = list(blocked_ids) + list(bound_blocked)
        non_blocked_ids = list(rxn_ids[blocked!=1])
        
        g = np.unique(fctable==1,axis=0)
        df = pd.DataFrame(data=g, columns= non_blocked_ids)
        coupled_sets=dict() #key = id of one reaction of the coupled set (first in alpha order); 
                            #value = list of coupled rxns ids
        for x in df.index:
            coupled = list(df.columns[df.loc[x]])
            coupled.sort()
            coupled_sets[coupled[0]] = coupled   
                
        self.coupled_rxns[prefix] = coupled_sets
        self.member_blocked[prefix]= blocked_ids
        total_rxns = len(rxn_ids)+len(bound_blocked)
        print("Flux coupling results for member %s stored:" % prefix)
        print("   Total reactions: %d" % total_rxns)
        print("   Fully coupled reaction sets: %d" % len(coupled_sets))
        print("   Blocked reactions: %d" % len(blocked_ids))      
        print("-")
        
    def set_cluster_reactions(self):
        
        # if FCA has been performed and results stored for all members, reactions for fva and clustering are reduced accordingly.
        # Otherwise, non-blocked reactions are obtained and ALL those reactions are used.       

        
        coupled_dicts =  list(self.coupled_rxns.values())
              
        if coupled_dicts.count(None) != 0: # at least one member without FCA results       
            print("Missing FCA results")
            print("Using non-blocked reactions only")
            self.get_non_blocked_reactions()  
            rxn2cluster = list(self.non_blocked)  
            #rxn2cluster = [r.id for r in self.cmodel.reactions]

        else:      
            accounted = []      
            coupled_rep = []
              
            # blocked reactions are not considered for fva and clustering   
            for prefix in self.member_blocked:
                blocked = self.member_blocked[prefix]
                accounted = accounted + blocked
              
            # Only one reaction from each fully coupled set is used for fva and clustering      
            for prefix in self.coupled_rxns:
                coupled_sets_dict = self.coupled_rxns[prefix]
                # rxn representative for coupled set:   
                coupled_rep  = coupled_rep + list(coupled_sets_dict.keys())
                coupled_sets = list(coupled_sets_dict.values())
                accounted = accounted + [rxn for cset in coupled_sets for rxn in cset]
        
            # Reactions for fva and clustering are those representing coupled sets and those not in any member (pool reactions)                   
            all_rids = [r.id for r in self.cmodel.reactions]
            missing = set(all_rids).difference(set(accounted))
        
            rxn2cluster = list(missing) + coupled_rep       
        
        rxn2cluster.sort()
        self.rxn2cluster =  rxn2cluster
        print("Total reactions considered for fva and clustering: %d" % len(self.rxn2cluster))
        
    def get_cluster_reaction_values(self, vector = 'qual_vector', thr=0.75, changing= True, convert=True):
        
        if self.clusters is None:
            raise RuntimeError("Missing clustering/qualitative FVA results!")
            
            

        
        #function to get representative qualitative values of a reaction in a cluster
        def get_rep_vals(x,thr):
            rep_val = None
            total = len(x)
            qual_vals, counts = np.unique(x, return_counts=True)
            rep = qual_vals[counts/total >= thr]
         
            if rep.size > 0:
                return rep[0]   #qualitative value present if at least thr of reactions in cluster  
            return None       
               
        
        #def changing_rep_vals(x):
        #    aux = x.unique()
        #    return aux.size > 1 
        
        #NJT to use qual_vector as well as bin_vector if required
        if vector =='qual_vector'and self.qual_vector_df is not None:
            z = self.qual_vector_df
        elif vector == 'bin_vector' and self.bin_vector_df is not None:
            z = self.bin_vector_df
        
        vector_df = z.astype('int32')
        
        cluster_ids = np.arange(1,self.k+1)
        cluster_dfs = [vector_df[self.clusters == c] for c in cluster_ids]
        aux = [ df.apply(get_rep_vals, thr=thr) for df in cluster_dfs]
        reps = pd.concat(aux, axis=1)
        reps.columns = ['c%d' % x for x in cluster_ids]

        reps = reps.astype('float')
        
        if changing: #report only reactions that have different qualitative values in at least two clusters
            changing_filter = reps.apply(lambda x: x.unique().size > 1, axis = 1)    
            reps = reps[changing_filter.values]
        
        if convert:
            cat_dict = {-3.0: '-', -2.0: '--',-1.0: '-0',1.0: '0+',0.0: '0',2.0: '++',3.0: '+',4.0: '-+',5.0: 'err',100.0: 'var'}
            reps = reps.replace(cat_dict)
        
        return reps
        
    @staticmethod
    def compare_clusters(clusters_df, cid_1, cid_2):
        
        #juice
        if isinstance(cid_1, int):
            cid_1 = 'c%d' % cid_1
        if isinstance(cid_2, int):
            cid_2 = 'c%d' % cid_2            
        
        df = clusters_df[[cid_1,cid_2]]
        changing_filter = df.apply(lambda x: x.unique().size > 1, axis = 1) 
        df = df[changing_filter.values]
        return df
        
    @staticmethod
    def plot_cluster_distribution(clusters_df, cmap='Accent',figsize=(10,5)):
        cl_columns = list(clusters_df)
        cat_percents_dict = dict()
        nreactions = clusters_df.shape[0]
        nan_rep = 100.0
        #1. Change nan to additional category 'variable'
        df2 = clusters_df.fillna(nan_rep)
        
        
        for c in cl_columns:
            vc = df2[c].value_counts()
            
            vc = 100 * vc/nreactions #to percentages
            cat_percents_dict[c] =  vc.to_dict()   
            
        cat_percents = pd.DataFrame.from_dict(cat_percents_dict, orient='index')
        cat_percents.fillna(0, inplace=True)
        
        cat_dict = {-3.0: '-',
                    -2.0: '--',
                    -1.0: '-0',
                     1.0: '0+',
                     0.0: '0',
                     2.0: '++',
                     3.0: '-+',
                     4.0: 'err',
                     5.0: '+',
                     100.0: 'var'}
        cat_percents.rename(columns = cat_dict, inplace=True)
        #plot
        ax = cat_percents.plot.barh(stacked=True, cmap=cmap,figsize=figsize)
        ax.legend(loc='center left',bbox_to_anchor=(1.0, 0.5),title='reaction category');
        ax.xaxis.set_major_formatter(mtick.PercentFormatter());
        ax.set_ylabel('clusters');
        ax.set_xlabel('reactions');

        return cat_percents 
    
    def clusterReactions(self, method, changing= True, **kwargs):
        cmodel = self.cmodel

        if self.qual_vector_df is None:
            print("No qualitative FVA values stored. Run qual_fva analysis first!")
            return

        print("Calculating jaccard distances between grid points...") 
        
        z = self.qual_vector_df.copy()
        self.changed_rxns = None
        if changing: #clustering considering changing reactions only:
            changed_rxns = self.qual_vector_df.max(axis=0) != self.qual_vector_df.min(axis=0)
            changed_rxns_ids = z.columns[changed_rxns]
            z = z[changed_rxns_ids]
            self.changed_rxns = changed_rxns_ids
            
        z = z.values
        z = z.T
        nrxns = z.shape[0]
        distance_metric = 'jaccard'
        dvector = distance.pdist(z,distance_metric)
        dmatrix = distance.squareform(dvector)
        
        print("Clustering %d reactions ..." % nrxns) 
        # Clustering methods
        # hierarchical + fcluster
        # dbscan (eps,min_samples)
        # optics
        # AffinityPropagation
        # SpectralClustering 
       
        if method == 'hierarchical':
            rk, rclusters = getHIERARCHICALclusters(dvector,**kwargs)
        elif method == 'dbscan':
            rk, rclusters = getDBSCANclusters(dmatrix,**kwargs)  
        elif method == 'optics':
            rk, rclusters = getOPTICSclusters(dmatrix,**kwargs)
        elif method == 'SpectralClustering':
            rk, rclusters = getSCclusters(dmatrix,**kwargs)          
        elif method == 'AffinityPropagation':
            rk, rclusters = getAPclusters(dmatrix, **kwargs)
            
        self.rk = rk
        self.rclusters = rclusters    
        print("Done!")    

    def quan_FCA(self, grid_x, grid_y, rxns_analysis):
        #Performs quantitative Flux Coupling Analysis on two reactions (rxns_analysis) and on points of a sub-grid defined by points grid_x, grid_y
        #returns: a dataframe with columns to plot afterwards
        #Columns: flux_rxns_analysis[0], flux_rxn_analysis[1], FVA (str: minimum or maximum), point (coordinates of point)

        feasible_points = self.points[self.feasible_points]
        analyze_points = []
        print('Quantitative Flux Coupling analysis \n Initializing grid...')
        def fraction_to_normalize(point_fractions, reaction):
            #from point_fraction computes which element of this array should be used for normalization
            #reaction: string reaction id
            fraction = ''
            for i, pre in enumerate(self.prefixes):
                if reaction.startswith(pre+'_'):
                    fraction = point_fractions[i]

            
            if fraction=='':
                print('No org detected, asumming community reaction')
                fraction =1
        
            return(fraction)
        #Match points defined by the user in grid_x, grid_y to specific points on the grid
        for y in grid_y:
            for x in grid_x:
                search_point = [x, y]
                distances = np.linalg.norm(feasible_points-search_point, axis=1)
                min_index = np.argmin(distances)
                analyze_points.append(min_index)
                print(f"the closest point to {search_point} is {feasible_points[min_index]}, at a distance of {distances[min_index]}")


        

        maxmin_data = []
        for this_point in analyze_points:
            model = copy.deepcopy(self)
            feasible_points = model.points[model.feasible_points]
            this_point_coords = feasible_points[this_point]
            print('Selected point'+str(this_point_coords))
            print('This point coords '+str(this_point_coords))
            this_point_frac = [this_point_coords[0], 1-this_point_coords[0]]
            print('This point frac '+str(this_point_frac))
            point = [this_point_coords[0]*this_point_coords[1], (1-this_point_coords[0])*this_point_coords[1]] #equivalent to old grid
            print('Old grid point '+str(point))

            #update bounds
            for i, member in enumerate(model.prefixes):
                mfrac = this_point_frac[i]
                mrxns = model.member_rxns[member]

                for rid in mrxns:
                    r = model.cmodel.reactions.get_by_id(rid)
                    old_bounds = r.bounds
                    r.bounds = (old_bounds[0]*mfrac, old_bounds[1]*mfrac)

            for ix, member_objectives in enumerate(model.objectives):
                new_bounds = (point[ix], point[ix])

                for rid in member_objectives.keys():
                    rxn = model.cmodel.reactions.get_by_id(rid)
                    rxn.bounds = new_bounds

            #try:
            #define limits reactions based on theoretical max-min defined from model
            rxn_ref_fva = flux_variability_analysis(model.cmodel, reaction_list = rxns_analysis[0])

            #define range reactions
            values_rxn_ref = np.linspace(rxn_ref_fva['minimum'][0], rxn_ref_fva['maximum'][0], num=50)
            values_rmax = []
            values_rmin = []

            with model.cmodel as cmodel:
                for val in values_rxn_ref:
                    rxn = cmodel.reactions.get_by_id(rxns_analysis[0])
                    rxn.bounds = (val,val)
                    #compute max min
                    fva = flux_variability_analysis(cmodel, reaction_list = rxns_analysis[1])
                    for i, el in enumerate(fva):
                        row_dict = dict()
                        row_dict[rxns_analysis[0]] = val/fraction_to_normalize(this_point_frac, rxns_analysis[0])
                        row_dict[rxns_analysis[1]] = fva[el][0]/fraction_to_normalize(this_point_frac, rxns_analysis[1])
                        row_dict['FVA'] = el
                        row_dict['point'] = str([round(this_point_coords[0],3), round(this_point_coords[1],3)])
                        maxmin_data.append(row_dict)

            #except:
            #    print('\n Issues with '+str(this_point_coords)+' unfeasible?')
        
        self.qFCA = pd.DataFrame(maxmin_data)
        

    def plot_qFCA(self, col_wrap=4):
        #Plots results computed by quan_FCA
        #input: maxmin_df (output of quan_FCA)
        #output: plot
        maxmin_df = self.qFCA
        sns.set(font_scale = 2)
        rxns_analysis = maxmin_df.columns[0:2]
        sns.set_style("whitegrid")

        g=sns.relplot(data = maxmin_df, x=rxns_analysis[0], y=rxns_analysis[1], col = 'point', hue='FVA', kind='line', col_wrap=4, lw=0)
        points = maxmin_df.point.unique()
        for i,ax in enumerate(g._axes):
            p = points[i]

            p_df = maxmin_df.loc[maxmin_df['point']==p]
            x = p_df.loc[p_df['FVA']=='maximum'][rxns_analysis[0]].to_numpy()

            y1 = p_df.loc[p_df['FVA']=='maximum']
            y1 = y1[rxns_analysis[1]].to_numpy()

            y2 = p_df.loc[p_df['FVA']=='minimum']
            y2 = y2[rxns_analysis[1]].to_numpy()

            ax.fill_between(x, y1,y2, color='none',hatch='//', edgecolor="k", linewidth=0.001)

