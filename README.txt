HOPF DELAY TOOLBOX

Cite: 
Synchronization in the Connectome: Metastable oscillatory modes emerge from interactions in the brain spacetime network
Joana Cabral, Francesca Castaldo, Jakub Vohryzek, Vladimir Litvak, Christian Bick, Renaud Lambiotte, Karl Friston, Morten L Kringelbach, Gustavo Deco

Contacts: 
joanacabral@med.uminho.pt
francesca.castaldo.20@ucl.ac.uk



 
FILE DESCRIPTION: total of 10 functions/scripts.

Part A- SIMULATE. Function 1) can be used to generate simulations and function 2) to extract spectral properties of simulated time series. In case one wants to run the whole parameter space running them on a cluster is highly reccomended. For downloading simulations, please contact the authors. 

Part B- ANALYSE and VISUALISE. To check/look at results/plots (DEMO results available)

Part C- DEMO Data. Structural Connectivity matrices to generate simulations; Model Spectral features (Metastability, Synchrony and Global Peak frequency).
MEG Power Spectra Density (PSD) for each subject and channel and the mean across channels for each subject. 
MEG Fitting matrices (Simulations PS for each pair of parameters, Squared Euclidean Distance matrices, Correlation Matrices between simulated PS and MEG PS). 
MEG HCP original data are not shared within the folder. For info and how do download them check the Manual: https://humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf

*************************************************************************************************************************************
Part A- SIMULATE.

1) hopf_delays_simu(f,K,MD,SynDelay,sig,varargin) *Gives Zsave as output (simulations) for each pair of the free parameters*

2) hopf_delays_feature_extraction(a,dt_save,MD,expK,C) *Loads simulations and extracts features (PSD, Metastability, Synchrony, Spectral Entropy [...])*
   2.1) bandpasshopf 
     2.1.1) convert_back_to_time 
*************************************************************************************************************************************

Part B- ANALYSE and VISUALISE.

1) network_parameter_space(a,MD,expK,C)
   *Load DEMO simu: Model_Spectral_Features.mat and create plots*

2) psd_meg_sensor_fit(a,dt_save,MD,expK,C)
  *Loads MEG Empirical (Sensor space), calculates and saves all simulated PSD, measures the model performance (fit by means of Squared Euclidean Distance and Correlation)*
   Output: MEG_sensor_PSD_Fitting.mat

3) meg_psd_model_screening(a, MD, expK, C) *This function loads simulations (MEG_sensor_PSD_Fitting.mat) and genrates model disparity plot (Squared Euclidean Distance) both for    each subject and averaged. 

4) MOM_analysis.mat 
   *Loads the baseline (no delay, optimal coupling in our ms), measures MOMs (5std above the baseline) for each band, measure MOMs Size, Duration, Occupancy; plots MOMs in all       brain areas for a given optimal point and bars for size, duration, occupancy*
        4.1) subplot_tight.mat *better than "subplot" - it helps in generating tight figures avoiding too much space between plots*

5) ModesInBands.mat
*Generates covariance matrices for selected points in the parameter space (specifically for No delays, weak coupling, Intermediate coupling - Optimal range, 
Strong Coupling, Long Delay) for each band, and how many frequency-specific eigenvalues are found for each of the points aforamentioned.*

***************************************************************************************************************************************

Part C- DEMO Data.

1) SC matrices 
File Name: SC_90aal_32HCP in Structural_Connectivity Folder
Citation: 

2) Simulation examples used to generate plots in the paper with a=-5; C=AAL90n32s; f=40Hz, MD=0:1:20; expK=-1:0.1:1.7; 
   Please download them here: https://liveuclac-my.sharepoint.com/:f:/g/personal/skgtfca_ucl_ac_uk/EjKFAcFpXC1FtJoCKVZoI1YBCB1gMZfZgl2SP83Tb9y9OA?e=1m7WF8 

For other simulations (i.e. different parcellation, different bifurcation parameter, different natural frequency) please write to francesca.castaldo.20@ucl.ac.uk

3) MEG DATA: https://humanconnectome.org/storage/app/media/documentation/s1200/HCP_S1200_Release_Reference_Manual.pdf


