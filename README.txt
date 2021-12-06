HOPF DELAYED MODEL TOOLBOX 

Cite: 
Synchronization in the Connectome: Metastable oscillatory modes emerge from interactions in the brain spacetime network
Joana Cabral, Francesca Castaldo, Jakub Vohryzek, Vladimir Litvak, Christian Bick, Renaud Lambiotte, Karl Friston, Morten L Kringelbach, Gustavo Deco

Contacts: 
joanacabral@med.uminho.pt
francesca.castaldo.20@ucl.ac.uk

 
FILE DESCRIPTION: total of 10 functions/scripts.

Part A- Run Simulations. Function 1 and 2 can be used to generate simulations and specrtral properties of simulated time series. 
However, this is not recommended unless you are considering this model for further investigation or to try different values of parameters. 

Part B- Model Analysis. To check/look at results/plots for the simulations available (DEMO simulations)

Part C- Data. Structural Connectivity matrices to generate simulations; Simulations used to generate the figures in the paper;  
MEG Power Spectra (PS) for each subject and channel and the mean across channels for each subject. MEG Fitting matrices (Simulations PS for each pair of parameters, Squared Euclidean Distance matrices, 
Correlation Matrices between simulated PS and MEG PS). MEG HCP original data are not shared within the folder, 
for info check Manual: https://humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf .


Part A- Run Simulations.

1) Hopf_delays_features.mat 
*Loads structural data (see SC matrices below for details of the example used in the paper), calls the next function (1.1) and measures features useful for the network dynamical analysis (Global Peak frequency, Metastability, Synchrony)*
  1.1) Hopf_Delays_Simu.mat  *Gives Zsave as output (simulations) for each pair of the free parameters*
  1.2) bandpasshopf 
       1.2.1)convert_back_to_time 

2) Hopf_delays_feature_extraction (useful if one wants to investigate additional features once simulations are available) 
*Loads simulations and extracts features (PSD, Metastability, Synchrony, Spectral Entropy [...])*


Part B- Model Analysis.

1) Network_Parameter_Space.mat 
*Load DEMO simu: Model_Spectral_Features.mat and create plots*

2) MEG_Fitting.mat 
*Loads MEG spectral features from HCP MEG 89 subjects (Power Spectrum) and simulations features (Power Spectrum) (load 'MEG_Fitting.mat'). Then, it measures and plots the model performance (Squared Euclidean Distance and Correlation)* 

3) MOM_analysis.mat 
*Loads the baseline (no delay, optimal coupling in our ms), measures MOMs in one optimal point (simu) (5std above the baseline) for each band, measure MOMs Size, Duration, Occupancy; plots MOMs in time in all brain areas for a given optimal point and bars for size, duration, occupancy; makes the video of MOMs in the optimal point over time*
3.1) subplot_tight.mat
*better than "subplot" - it helps in generating tight figures avoiding too much space between plots*

4) ModesInBands.mat
*Generates covariance matrices for selected points in the parameter space (specifically for No delays, weak coupling, Intermediate coupling - Optimal range, 
Strong Coupling, Long Delay) for each band, and how many frequency-specific eigenvalues are found for each of the points aforamentioned.*


Part C- Data.

1) SC matrices 
File Name: SC_90aal_32HCP in Structural_Connectivity Folder
Citation:  

2) Simulation examples used to generate plots in the paper. Please download them here: 
https://liveuclac-my.sharepoint.com/:f:/g/personal/skgtfca_ucl_ac_uk/Etzx5-w3LFhMgPhJGXBErlgBwxCK2VYyTOHgUNk0pivHtw?e=kgmcgn  



