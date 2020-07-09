# Biroegion_Methods
Comparison of analytical methods for bioregionalisation

This repository contains the R Code used for the Methods in Ecology and Evolution paper, Determining Marine Bioregions: A comparison of quantitative approaches by Hill et al. 2020.

It consists of two main folders, one for the generation and analysis of simulation data and the other for the analysis of the Kerguelen Plateau demersal fish data.
The contents of the folders are described below.

Simulation:

The ‘Simulation_Setup’ folder contains the files used to produce the simulation data.
- simulation_env_070518.R: Code to generate environmental covariates across a simulated region.
- sim_env_070518.RData: The environmental covariates across the simulated region generated from the above code.
- simulate_communities_final.R: Code to generate simulated communities comprising the bioregions.
- Many_covars_sim_fin.RData: The simulated bioregion outputs.

The remaining files contain code to run the various bioregionalisation methods and interpret their outputs.
- Simulation_Run_Models.R: Runs the bioregionalisation analyses.
- Simulation_Additional_Funcs.R: Contains extra functions used to run models and interpret results.
- Plot_pred_clusters.R: Code to plot the predicted distribution of bioregions.
- Plot_Group_Contents.R: Code to generate and plot the species composition of bioregions.
- Plot_environment.R: Code to generate and plot the environmental characteristics of bioregions.


Files for the analyses of the Kerguelen Plateau fish dataset are contained in the folder 'KP_Fish':
- KPFish_Run_Models.R: Runs the bioregionalisation analyses.
- Additional_Funcs.R: Contains extra functions used to run models and interpret results.
- KPFish_Plot_pred_clusters.R: Code to plot the predicted distribution of bioregions.
- KPFish_Plot_Group_Contents.R: Code to generate and plot the species composition of bioregions.
- KPFish_Plot_environment.R: Code to generate and plot the environmental characteristics of bioregions.

KP Fish data and associated environmental predictor variables needed to run the above analysis code can be found here:
http://dx.doi.org/doi:10.26179/5f0528de8c1d2

Rasters of the environmental variables and the prediction space needed to generate spatial predictions of the distribution of bioregions can be found here:
http://dx.doi.org/doi:10.26179/5f055cd217aa8
