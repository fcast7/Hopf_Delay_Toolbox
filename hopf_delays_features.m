function hopf_delays_features

%%%%%%%%%%%%%%%%%%%%%%%
%  RUN SIMULATIONS 
%  This is the function to extract Spectral features from simulated data
%  Zsave
%  
%
%  Input:   Structural connectivity (fiber length and distances)
%
%  Output:  Synchrony, Metastability, Global Peak Frequency, Mean Frequency,
%           Global Spectral Entropy 
% 
%  Nested functions (add to path): 
%  - Hopf_Delays_Simu.mat
%  - bandpasshopf.mat 
%  - convert_back_to_time.mat
%
%  Joana Cabral 2017 joana.cabral@psych.ox.ac.uk
%  Francesca Castaldo 2021 francesca.castaldo.20@ucl.ac.uk
%  
%%%%%%%%%%%%%%%%%%%%%%%%

% Define the Structural Network (all fibers) 
addpath('Hopf_Delay_Toolbox');
load SC_90aal_32HCP.mat mat
N=size(mat,1);
C=mat/mean(mat(ones(N)-eye(N)>0)); %% Such that the mean of all non-diagonal elements is 1.

%Reduce the length to speed-up code. 
%Define the Structural Network (only number of fibers tracts detected between 2ROIs>10, so when we cut the connections
%with less fibers, we are mostly removing the longest fibers and not the shortest) 
% load SC_90aal_32HCP.mat mat
% red_mat=mat;
% red_mat(red_mat<10)=0;
% N=size(red_mat,1);
% C=red_mat/mean(red_mat(ones(N)-eye(N)>0));
% Such that the mean of all non-diagonal elements is 1.

% Mean Connection distance
% v=mean(D(C>0))/MD 


% Distance between areas
load SC_90aal_32HCP.mat mat_D
D=mat_D;
D=D/1000; % Distance matrix in meters

% Define simulation Parameters
tmax=50; % time of simulation in seconds
t_prev=5; % in seconds
dt_save=2e-3; % max resolution to save

f=40; %Intrinsic frequency in Hz
sig=1e-3; % noise std 1e-3

% First Free Parameter: Mean Delay 
MD=0:1:20; % Reange of Mean Delay in ms or a fixed MD
MD=MD.*1e-3; % Mean Delay in seconds 
SynDelay=0; % Constant Synaptic Delay (Do not consider it if the model still need to be tested)

% Second Free parameter: coupling strength 
expK=-1:0.1:1.7; %The max limit of this range depens also on the dt - issue: amplitude death 
K=10.^(expK);


Sync               = zeros(length(K),length(MD));
Meta               = zeros(length(K),length(MD));
PeakFMean          = zeros(length(K),length(MD));
PeakFGlobal        = zeros(length(K),length(MD));
SpecEntropy_Global = zeros(length(K),length(MD));


fbins=1000;
freqZ = (0:fbins-1)/(dt_save*fbins);

for g=1:length(K)
    for d=1:length(MD)
        k=K(g);
        md=MD(d);
        
        [Z] = hopf_delays_simu(f,k,md,SynDelay,sig,C,D,tmax,t_prev,dt_save);
                
        % Detect the peak Frequency in the Fourier Transform of all areas        
        Fourier_Complex = fft(Z,fbins,2); %% Fourier of Z (complex) in 2nd dimension
                         
        Fourier_Global=abs(mean(Fourier_Complex)).^2;
        [~, Imax]=max(Fourier_Global);
        PeakFGlobal(g,d)=freqZ(Imax);
        
         % Evaluate Order Stability around the Global peak frequency
         % Band-pass filter the signals Z around the peak frequency of the
         % ensemble 
        for n=1:90
            Zb(n,:)=bandpasshopf(Z(n,:),[max(0.1,freqZ(Imax)-1) freqZ(Imax)+1],1/dt_save);
            Zb(n,:)=angle(hilbert(Zb(n,:)));
        end
        OP =abs(mean(exp(1i*(Zb)),1));
        
        Sync(g,d)=mean(OP); % measure of global synchronization 
        Meta(g,d)=std(OP); %  how much R (the order parameter) fluctuates in time      

             
        % Probability distribution (total sum is 1)
        PSD_G=Fourier_Global/sum(Fourier_Global);
        
        %Spectral Entropy method 1
        PSD=PSD_G/sum(PSD_G+1e-12);
        logPSD = log2(PSD_G + 1e-12);
        SpecEntropy_Global(g,d) = -sum(PSD_G.*logPSD)/log2(length(PSD_G)); %Normalized spectral entropy 
        
     
        % Average of Node PSD
        Fourier_Mean=mean(abs(Fourier_Complex).^2);
        [~, Imax]=max(Fourier_Mean);
        PeakFMean(g,d)=freqZ(Imax);
        
        % Probability distribution (total sum is 1)
        PSD_M=Fourier_Mean/sum(Fourier_Mean);

        K_label=num2str(log10(K));
        ind_p=find(K_label=='.');
        if numel(ind_p)
        K_label(ind_p)='p';
        end
%         save(['AAL_HCP' num2str(MD*1e3) 'MD_' K_label 'K_' num2str(f) 'Hz' ],'MD','K','Sync','Meta','PeakFGlobal','PeakFMean', 'PSD_G','PSD_M')
        save(['Hopf_delays_features' num2str(f) 'Hz' ],'MD','K','Sync','Meta','PeakFGlobal','PeakFMean','SpecEntropy_Global')

    end    
end

