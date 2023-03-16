from sklearn.cluster import DBSCAN, OPTICS,SpectralClustering, AffinityPropagation

from scipy.cluster.hierarchy import fcluster
from scipy.cluster import hierarchy

from itertools import compress

import cobra.io
import pickle
import numpy as np
import pandas as pd

#  Funciones para correr en paralelo en el cluster y luego cargar resultados

def store_miniEco_inputs(eco, file_dir, name_scenario,name_community_model):
    
    
    #cobra model for community
    
    cmodel_filename = '%s/%s_sce_%s.json' % (file_dir,name_community_model,name_scenario)  
    cobra.io.save_json_model(eco.cmodel, cmodel_filename)        
    
    #relevant metadata for fva and qual transform
    if eco.clusters is not None:
        data_dict = {'points': eco.points,
                 'objectives':eco.objectives,
                 'prefixes':eco.prefixes,
                 'rxn2cluster':eco.rxn2cluster,
                 'member_rxns':eco.member_rxns,
                 'feasible_points':eco.feasible_points,
                 'clusters': eco.clusters}
    else:
        data_dict = {'points': eco.points,
                 'objectives':eco.objectives,
                 'prefixes':eco.prefixes,
                 'rxn2cluster':eco.rxn2cluster,
                 'member_rxns':eco.member_rxns,
                 'feasible_points':eco.feasible_points}

    
    pickle_filename = '%s/%s_sce_%s.p' % (file_dir,name_community_model,name_scenario) 
    pickle.dump(data_dict, open( pickle_filename, "wb" ) )
        
    print("Reduced data for model %s and scenario %s stored." % (name_community_model, name_scenario))
    print(" metadata output: %s" % pickle_filename)                          
    print(" cobra model output: %s" % cmodel_filename)   

def check_batch_files(eco,file_dir,file_prefix):
    
    # Check if points and reactions2cluster are the same as those stored in pickle_file
    pickle_filename = '%s/%s.p' % (file_dir,file_prefix) 
    data_dict = pickle.load( open( pickle_filename, "rb" ) )
    
    check1 = eco.rxn2cluster != data_dict['rxn2cluster']
    if check1:
        aux1 = sorted(eco.rxn2cluster)
        aux2 = sorted(data_dict['rxn2cluster'])
        if aux1 == aux2:
            print("Warning rxns stored in different order") 
            check1 = False
    
    check2 = not np.array_equal(eco.points,data_dict['points'])
    if check1 or check2 :
        print("Different data used to generate qualitative Matrix")
        return False        
    return True     

def get_feasible_points_from_batch_files(eco, file_dir,name_scenario,name_community_model, batch_size):
    
    file_prefix = '%s_sce_%s' % (name_community_model,name_scenario) 
     #correct_files = check_batch_files(eco,file_dir,file_prefix)
    #if not correct_files:
        #raise RuntimeError("Wrong input files")

    #Build feasible point boolean list from batch files 
    npoints = len(eco.points)
    full_batches = npoints // batch_size
    extra_batch = npoints % batch_size         
    batches = full_batches + (extra_batch > 0) *1
    batches_first = [x*batch_size for x in range(batches)]
    batch_files = ["%s/%s_fpoints_%d_%d.npy" % (file_dir, file_prefix, first ,(first + batch_size)) 
                   for first in batches_first]
    if extra_batch > 0:
        batch_files[-1] = "%s/%s_fpoints_%d_%d.npy" % (file_dir, file_prefix, batches_first[-1] ,npoints) 
    batch_vectors = [np.load(file) for file in batch_files]    
    feasible_points = list(np.concatenate(batch_vectors))

    return feasible_points
    
def get_qual_fva_from_batch_files(eco, file_dir,name_scenario,name_community_model, batch_size):
   
    file_prefix = '%s_sce_%s' % (name_community_model,name_scenario) 
    #correct_files = check_batch_files(eco,file_dir,file_prefix)
          
    #if not correct_files:
    #    raise RuntimeError("Wrong input files")        
   
    #batches indexes
    npoints = len(eco.points)
    full_batches = npoints // batch_size
    extra_batch = npoints % batch_size
    batches = full_batches + (extra_batch > 0) *1
    batches_first = [x*batch_size for x in range(batches)]

    if eco.feasible_points is None: 
        # we assume FVA was run on full grid including unfeasible community steady-state points
        filename_flag = 'qFVA' 
        df_index = range(len(eco.points))
        valid_batches = [True for x in range(len(batch_files))]
        
    else:
        filename_flag = 'qFVA_feasibleOnly'
        df_index = np.where(eco.feasible_points)[0]  
        #check for batches with feasible points only:
        valid_batches = [sum(eco.feasible_points[x:min(x+batch_size,npoints)])>0 for x in batches_first]  
        
    #qualitative FVA batch files
    batch_files = ["%s/%s_%s_%d_%d.npy" % (file_dir, file_prefix, filename_flag,first ,(first + batch_size)) for first in batches_first]
    if extra_batch > 0:
        batch_files[-1] = "%s/%s_%s_%d_%d.npy" % (file_dir, file_prefix, filename_flag, batches_first[-1] ,npoints)
    batch_files = list(compress(batch_files,valid_batches))
   	
    batch_arrays = [np.load(file) for file in batch_files]    
    full_array = np.concatenate(batch_arrays)

    qual_vector_df = pd.DataFrame(data = full_array, columns = eco.rxn2cluster, index = df_index)
    #qual_vector_df = pd.DataFrame(data = full_array, columns = data_dict['rxn2cluster'], index = df_index)    
     
    return qual_vector_df   


def get_fva_from_batch_files(eco, file_dir,name_scenario,name_community_model, batch_size):
   
    file_prefix = '%s_sce_%s' % (name_community_model,name_scenario) 
    #correct_files = check_batch_files(eco,file_dir,file_prefix)
          
    #if not correct_files:
    #    raise RuntimeError("Wrong input files")        
   
    #batches indexes
    npoints = len(eco.points)
    full_batches = npoints // batch_size
    extra_batch = npoints % batch_size
    batches = full_batches + (extra_batch > 0) *1
    batches_first = [x*batch_size for x in range(batches)]

    if eco.feasible_points is None: 
        # we assume FVA was run on full grid including unfeasible community steady-state points
        filename_flag = 'full_FVA'
        df_index = range(len(eco.points))
        valid_batches = [True for x in range(len(batch_files))]
        
    else:
        filename_flag = 'full_FVA_feasibleOnly'
        df_index = np.where(eco.feasible_points)[0]  
        #check for batches with feasible points only:
        valid_batches = [sum(eco.feasible_points[x:min(x+batch_size,npoints)])>0 for x in batches_first]  
        
    #qualitative FVA batch files
    batch_files = ["%s/%s_%s_%d_%d.npy" % (file_dir, file_prefix, filename_flag,first ,(first + batch_size)) for first in batches_first]
    if extra_batch > 0:
        batch_files[-1] = "%s/%s_%s_%d_%d.npy" % (file_dir, file_prefix, filename_flag, batches_first[-1] ,npoints)
    batch_files = list(compress(batch_files,valid_batches))
   	
    batch_arrays = [np.load(file) for file in batch_files]    
    full_array = np.concatenate(batch_arrays)

    #fva = pd.DataFrame(data = full_array, columns = eco.rxn2cluster, index = df_index)
    fva = full_array 
     
    return fva 



# Funciones para clusterizar puntos

def getHIERARCHICALclusters(dvector,k=20,lmethod='ward',criterion='maxclust', **kwards):
    row_linkage = hierarchy.linkage(dvector, method=lmethod)
    clusters = fcluster(row_linkage, k, criterion=criterion)
    return k, clusters

def getDBSCANclusters(dmatrix,eps = 0.05, min_samples = 5, **kwards): 
    
    # eps : 
    # The maximum distance between two samples for one to be considered as in the neighborhood of the other. 
    # This is not a maximum bound on the distances of points within a cluster. This is the most important 
    # DBSCAN parameter to choose appropriately for your data set and distance function.
 
    # min_samples:
    # The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. 
    # This includes the point itself.
    
    dbscan = DBSCAN(eps=eps, min_samples = min_samples,metric='precomputed',**kwards)
    clusters = dbscan.fit_predict(dmatrix) 
    
    #output incluye outliers, valor -1 en vector clusters
    k = len(np.unique(clusters))
    
    cluster_ids = list(np.arange(k)+1)
    cluster_db_ids = list(np.unique(clusters))
    cluster_dict = dict(zip(cluster_db_ids, cluster_ids))
    print("outliers: cluster %d" % cluster_dict[-1])
    clusters = [cluster_dict[x] for x in clusters]
    
    
    return k,clusters
    
def getOPTICSclusters(dmatrix, max_eps=0.05, min_samples = 5, **kwards):
    optics = OPTICS(max_eps=max_eps, min_samples = min_samples,metric='precomputed')
    clusters = optics.fit_predict(dmatrix)
    k = len(np.unique(clusters))

def getSCclusters(dmatrix,assign_labels= "discretize",random_state=0,k=20, delta=0.2, **kwards):

    
    #transformación de matriz de distancia a matriz de similitud. Vía aplicación de Gaussian (RBF, heat) kernel:
    sim_matrix = np.exp(- dmatrix ** 2 / (2. * delta ** 2))
    
    sc = SpectralClustering(n_clusters=k, assign_labels=assign_labels, 
                            random_state=random_state,affinity='precomputed', **kwards)
    sc_results = sc.fit_predict(sim_matrix)
    #clusters = sc_results.labels_
    clusters = sc_results
    return k, clusters   

def getAPclusters(dmatrix, delta=0.2, **kwards):
    af = AffinityPropagation(**kwards)
    #sim_matrix = np.exp(- dmatrix ** 2 / (2. * delta ** 2))
    af_matrix = dmatrix
    af_results = af.fit_predict(af_matrix)
    #clusters = af_results.labels_
    clusters = list(af_results)
    k = len(np.unique(clusters))
    return k, clusters
