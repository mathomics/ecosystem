# ecosystem
Supplementary material _"Unveiling abundance-dependent metabolic phenotypes of microbial communities"_, consisting on an framework that allows to explore the abundance-growth space, defined by both the composition of a community as well as its growth rate. By performing a cualitative characterization of each point of these space, this method allows to characterize metabolic phenotypes in a community without the imposition of a biomass maximization objective function. 

This framework is implemented in Python and uses the individual genome-scale models as an input to construct a community model in an `Ecosystem` object as described in the available tutorials on this repository. 

## Install instructions
Currently the raw code can be directly downloaded to be used from this repository
- `full_ecosystem.py`: class of the ecosystem object
- `multiple_ecosystem.py`: class of the ecosystem object optimized for multiple analysis (no model is loaded for each scenario for faster computation, only grid and results computed using `full_ecosystem.py`)

## Basic usage
A basic set of command that allow to perform these analysis is described bellow, for a more detailed description see `sup_toy_example.ipynb` which includes additional explanation and plots.

### Initiation and simulation conditions
To generate an 'Ecosystem' object for a community of two organisms with models `model_1` and `model_2`, which contains a constrained-based model of the community:
```
eco2 = Ecosystem(models=[model_1, model_2], prefixes=['org1', 'org2], community_name = 'my first community', community_id = 'com_2', solver=gurobi)
```
To constrain inputs of nutrients (in this case metabolite `A_e`) to the common space shared by all members of the community:
```
eco2.set_pool_bounds({'A_e':(-lb,1000)}, bioCons=-bc)
```
This function also accepts an additional parameter where biological constraints (bioCons) can be set up for the organisms of the community.

To generate a grid of 10x10 for analysis:
```
eco2.set_cluster_reactions()
eco2.build_grid(numPoints = 10, expand = True, drop_zero=True)
eco2.get_member_reactions()
eco2.get_points_distribution()
```
### Grid analysis and clustering
To analyze the grid different tests can be performed (feasibility, qualitative Flux Variability analysis (FVA)), as follows:
``` 
eco2.analyze_grid(analysis = 'feasibility', update_bounds=True) 
```
checks the feasibility of each point of the pre-computed grid of abundance and community growth rates
or 
```
eco2.analyze_grid(analysis = 'qual_fva', update_bounds=True)
```
which performs qualitative FVA on each point of the grid. The output of this analysis can be analyzed by clustering a defined number of clusters (`num_clusters`):

```
eco2.clusterPoints('hierarchical', k = num_clusters)
```
which can be plot using:
```
eco2.plot_2D_slice(prefixes=[], fixed_values=[], parent_cmap='tab20c',s=70, figsize=(12,12), 
                         to_plot = 'cluster', show_edge=False,frac_prefix= None,
                            xlabel = '$f_{org1}$',
                            ylabel ='$Biomass_{community}$', saveFile='clusters_toy')
```

or analyzed on a table format by using 

```
eco2.get_cluster_reaction_values(thr=0.8, changing= True)
```
 and
```
eco2.compare_clusters(df, 'c3','c4')
```
to compare two clusters (in this case `c3` and `c4`.

### Quantitative Flux Coupling Analysis
A quantitative version of Flux Coupling Analysis (FCA) can be performed on different points of the grid defined by their coordinates, and for two reactions of choice to study how distribution of fluxes changes on different points of the grid, as follows:
```
eco2.quan_FCA(grid_x, grid_y, rxns_analysis)
plot_qFCA(col_wrap=4)
```


## Supplementary files
- S1: Model curation for _A. ferrooxidans_ and _S. thermosulfidooxidans_
- S2: Models _A. ferrooxidans_ Wenelen (iML510), _S. thermosulfidooxidans_ Cutipay (iSM517) and the community model generated for the bioleaching community in the base case scenario, in both sbml and xls format.

## Tutorials and case examples
- `sup_toy_example.ipynb` includes an analysis that can be performed locally and its associated plots
- `sup_mat_multiple_bioleaching.ipynb` includes analysis of multiple scenarios simultaneously from previously computed results stored on the folder result_files. It requires both models included in Supplementary File 2,  `simulations_pareto4.csv` to filter the different simulated scenarios and the files available on the folder `result_files`

## How to cite
- Jiménez, NE., Acuña V., Cortés, MP., Eveillard, D., Maass A., (2023) Unveiling abundance-dependent metabolic phenotypes of microbial communities
