# ecosystem
Supplementary material _"The abundance-growth space determines relevant metabolic phenotypes of a community"_, which consists on an framework that allows to explore the abundance-growth space, defined by both the composition of a community as well as its growth rate. By performing a cualitative characterization of each point of these space, this method allows to characterize metabolic interactions in a community without the imposition of a biomass maximization objective function. 

This framework is implemented in Python and uses the individual genome-scale models as an input to construct a community model in an `Ecosystem` object as described in the available tutorials on this repository

## Install instructions
Currently the raw code can be directly downloaded to be used from this repository
- `full_ecosystem.py`: class of the ecosystem object
- `multiple_ecosystem.py`: class of the ecosystem object optimized for multiple analysis (no model is loaded for each scenario for faster computation, only grid and results computed using `full_ecosystem.py`)

## Supplementary files
- S1: Model curation for _A. ferrooxidans_ and _S. thermosulfidooxidans_
- S2: Models _A. ferrooxidans_ Wenelen (iML510) and _S. thermosulfidooxidans_ Cutipay (iSM517)

## Tutorials and case examples
- `sup_mat_multiple_bioleaching.ipynb` includes analysis of multiple scenarios simultaneously from previously computed results stored on the folder result_files. It required both models included in Supplementary File 2,  `simulations_pareto4.csv` to filter the different simulated scenarios and the files available on the folder `result_files`
- `sup_toy_example.ipynb` includes an analysis that can be performed locally and its associated plots

## How to cite
- Jiménez, NE., Acuña V., Cortés, MP., Eveillard, D., Maass A., (2022) The abundance-growth space determines relevant metabolic states of a community
