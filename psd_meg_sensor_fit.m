function psd_meg_sensor_fit(a,dt_save,MD,expK,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function to explore a range of parameters in the Whole-Brain Hopf Model
%  and compare with MEG Empirical data
%
%  Francesca Castaldo (francesca.castaldo.20@ucl.ac.uk) 2022
%  
% Input:
%
%  a: bifurcation parameter
%  tmax: time of the simulations
%  dt_save: simulation resolution
%  MD: range of mean delays in ms
%  expK: range of Coupling
%  C: structural connectivity
%
%  Example: psd_meg_sensor_fit(-5,2e-3,0:1:20,-1:0.1:1.7,AAL90n32s)
%
%      1) Define parameters
%      2) Loads MEG Empirical (Sensor space)
%      3) Calculates and saves all simulated PSD
%      4) Measures the model performance (Squared Euclidean Distance and Correlation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE the current folder to save the output

if a==-5
    if C==AAL90n32s
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg5')
        disp(['Running for' num2str(a) '90AAL32'])
    elseif C==AAL90n985s
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\985AAL\a_neg5')
        disp(['Running for' num2str(a) '90AAL985'])
    elseif C==SHEAF200n32s
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[SHEAFER]\Simulations')
        disp(['Running for' num2str(a) '200SHEAF32'])
    end

elseif a==-0.2
    if C==AAL90n32s
    cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg02')
    disp(['Running for' num2str(a) '90AAL32'])
    end
    
elseif a==-0.05
    if C==AAL90n32s
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg005\Simulations')
        disp(['Running for' num2str(a) '90AAL32'])
    end
end


% 1) Define simulation Parameters
K=10.^(expK);
MD=MD.*1e-3; % Mean Delay in seconds
fbins=1000;
freqZ = (0:fbins-1)/(dt_save*fbins);
% ind5Hz=find(freqZ>=5,1);
% ind25Hz=find(freqZ>=25);
% ind_Low_F=find(freqZ==3); %in case we want to investigate only frequency above this value
ind150Hz=find(freqZ==150); %range of interest


PSD_Simu_Global      = zeros(length(K),length(MD),ind150Hz); %whole signal
% PSD_SimuSensor    = zeros(length(K),length(MD),41); %5-25Hz signal
Error_MEG_PSD       = zeros(length(K),length(MD));
Fit_MEG_PSD         = zeros(length(K),length(MD));
Pdist_Error_MEG_PSD = zeros(length(K),length(MD));
Error_MEG_PSD_Sub   = zeros(length(K),length(MD), 89);
Corr_Fit_MEG_PSD    = zeros(length(K),length(MD));
% fit_MEG_envFC       = zeros(length(K),length(MD));
% corrSC              = zeros(length(K),length(MD));

% 2) Load Empirical Data
% Load Mean Across Channels for each subject (MEG_MeanPSD_Planar_89) and
% mean over subjects (Mean_PSD_Planar_89). Data can be found at the HCP
% website. 

load('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\MEG_Data\MEG_Sensor\Results\MEG_MeanPSD_Planar_89.mat')
load('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\MEG_Data\MEG_Sensor\Results\Mean_PSD_Planar_89.mat')


% 3) Calculate and save the model PSD in a file called MEG_PSD_Fitting (see
% below)

for g=1:length(K)
    for d=1:length(MD)
        k=K(g);
        md=MD(d);

        disp(['Now K=' num2str(k) ', mean Delay = ' num2str(md*1e3) 'ms'])

        K_label=num2str(log10(k));
        ind_p=find(K_label=='.');
        if numel(ind_p)
            K_label(ind_p)='p';
        end

        load(['a_Remote_K1E' K_label '_MD_' num2str(md*1e3) 'a' num2str(a)],'Zsave') %for cluster new simulations


        Fourier_Simu_Global= fft(mean(Zsave),fbins); %% Fourier of Z (complex)
        Simu_PSD_MEAN=abs(Fourier_Simu_Global(:,1:ind150Hz)).^2;
        Simu_PSD_MEAN=Simu_PSD_MEAN/sum(Simu_PSD_MEAN);

        %           Simu_PSD_MEAN=squeeze(PSD_Simu_Global(g,d,:))'; Uncomment this
        %           line and comment the lines up if you have PSD_Simu_Global
        %           already saved

        Error_MEG_PSD(g,d)=pdist([cumsum(Mean_PSD_Planar_89); cumsum(Simu_PSD_MEAN)],'squaredeuclidean');
        Fit_MEG_PSD(g,d)= corr2(cumsum(Mean_PSD_Planar_89), cumsum(Simu_PSD_MEAN));
        Pdist_Error_MEG_PSD(g,d)=pdist([Mean_PSD_Planar_89;Simu_PSD_MEAN],'squaredeuclidean');
        Corr_Fit_MEG_PSD(g,d)= corr2(Mean_PSD_Planar_89,Simu_PSD_MEAN);


        for s=1:89
            Subj_PSD_Planar=MEG_MeanPSD_Planar_89(s,:)/sum(MEG_MeanPSD_Planar_89(s,:));
            Error_MEG_PSD_Sub(g,d,s) = pdist([cumsum(Subj_PSD_Planar); cumsum(Simu_PSD_MEAN)]);%,'squaredeuclidean');
            Error_MEG_PSD_Sub_squared(g,d,s) = pdist([cumsum(Subj_PSD_Planar); cumsum(Simu_PSD_MEAN)],'squaredeuclidean');
            Corr_MEG_PSD_Sub(g,d,s)=corr2(cumsum(Subj_PSD_Planar), cumsum(Simu_PSD_MEAN));
        end
        PSD_Simu_Global(g,d,:) = Simu_PSD_MEAN;

    end
end

save('MEG_sensor_PSD_Fitting','MD','expK','Error_MEG_PSD','PSD_Simu_Global','Error_MEG_PSD_Sub','Corr_Fit_MEG_PSD','Fit_MEG_PSD','Pdist_Error_MEG_PSD','Corr_MEG_PSD_Sub','Error_MEG_PSD_Sub_squared')

