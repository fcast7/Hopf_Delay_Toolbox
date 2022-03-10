function hopf_delays_feature_extraction(a,dt_save,MD,expK,C) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Whole-brain Hopf Bifurcation Model with time delay%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%      This function loads simulated time series (Zsave - from hopf_delays_simu) and measures 
%      Spectral and global dynamical features on simulated data 
%
%  - Input: Define Simulation Parameter 
%  
%     
%     a:       bifurcation parameter 
%     dt_save: resolution or downsampling 
%     MD:      range of delays in ms 
%     expK:    range of couplings (expK)
%     C:       structural connectome, nodes_parcellation_nsubj (i.e. AAL90n32s) 
%
%   
%  - Output: Spectral features of simulated time series in the Whole-Brain Hopf Model
%             to save as Model_Spectral_Features.mat (uncomment fo measuring also Spectral
%             Entropy and Lempel-Ziv Complexity)
%
%  Examples: 
%  32AAL, -5 simu: Hopf_delays_feature_extraction(-5,2e-3,0:1:20,-1:0.1:1.7,AAL90n32s)
%  985AAL, -5 simu: Hopf_delays_feature_extraction(-5,2e-3,0:1:20,-1:0.1:1.5,AAL90n985s)
%  
%  Demo:
%  Hopf_delays_feature_extraction(-5,2e-3,0:1:5,-1:0.1:1.7,AAL90n32s)
%  
%  Code by Francesca Castaldo 2022 francesca.castaldo.20@ucl.ac.uk
%  Adapted from Joana Cabral 2017 joana.cabral@psych.ox.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Define your own path for accessing simulations
if a==-5
    if C=='AAL90n32s'
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg5')
        disp(['Running for' num2str(a) '90AAL32'])
    elseif C=='AAL90n985s'
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\985AAL\a_neg5')
        disp(['Running for' num2str(a) '90AAL985'])
    elseif C=='SHEAF200n32s'
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[SHEAFER]\Simulations')
        disp(['Running for' num2str(a) '200SHEAF32'])
    end

elseif a==-0.2
    if C=='AAL90n32s'
    cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg02')
    disp(['Running for' num2str(a) '90AAL32'])
    end
    
elseif a==-0.05
    if C=='AAL90n32s'
        cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg005\Simulations')
        disp(['Running for' num2str(a) '90AAL32'])
    end
end

% Define simulation Parameters Scaling
MD=MD.*1e-3; % Mean Delay in seconds
K=10.^(expK);
fbins=1000;
freqZ = (0:fbins-1)/(dt_save*fbins);

% Inizialize 
Sync                 = zeros(length(K),length(MD));
Meta                 = zeros(length(K),length(MD));
PeakFMean            = zeros(length(K),length(MD));
PeakFGlobal          = zeros(length(K),length(MD));
% SpecEntropy_Global = zeros(length(K), length(MD));
% lzcomplexity       = zeros(length(K), length(MD));

%%  2) Generate spectral features given simulated time series

% figure('color', 'w') %only if one wants to see the simulations in the whole parameter space

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
       
        % load simulations 
        load(['a_Remote_K1E' K_label '_MD_' num2str(md*1e3) 'a' num2str(a)],'Zsave')


        % Detect the peak Frequency in the Fourier Transform of all areas
        Fourier_Complex = fft(Zsave,fbins,2); %% Fourier of Z (complex) in 2nd dimension
       
        Fourier_Global=abs(mean(Fourier_Complex)).^2;
        [~, Imax]=max(Fourier_Global);
        PeakFGlobal(g,d)=freqZ(Imax);
        
        % Evaluate Order Stability around the Global peak frequency
        % Band-pass filter the signals Z around the peak frequency of the
        % ensemble
       
        for n=1:size(Zsave,1)
            Zb(n,:)=bandpasshopf(Zsave(n,:),[max(0.1,freqZ(Imax)-1) freqZ(Imax)+1],1/dt_save);
            Zb(n,:)=angle(hilbert(Zb(n,:)));
            
        end
        OP =abs(mean(exp(1i*(Zb)),1));
        
        Sync(g,d)=mean(OP); % measure of global synchronization (mean of OP)
        Meta(g,d)=std(OP); %  how much OP fluctuates in time (std of OP)
        
        % Probability distribution (total sum is 1)
%         PSD_G=Fourier_Global/sum(Fourier_Global);
        
        
        % Average of Node PSD
        Fourier_Mean=mean(abs(Fourier_Complex).^2);
        [~, Imax]=max(Fourier_Mean);
        PeakFMean(g,d)=freqZ(Imax);
        
        % Probability distribution (total sum is 1)
%         PSD_M=Fourier_Mean/sum(Fourier_Mean);
        
        %Spectral Entropy method 1
%         PSD=PSD_G/sum(PSD_G+1e-12);
%         logPSD = log2(PSD_G + 1e-12);
%         SpecEntropy_Global(g,d) = -sum(PSD_G.*logPSD)/log2(length(PSD_G)); %Normalized spectral entropy
        
        %         %Spectral Entropy Method 2 + lempel ziv complexity
        %         for n=1:90
        %           SE = pentropy(real(Zsave(n,:)), 1/dt_save, 'Instantaneous', false);
        %           LZC= LZComp_median(Zsave(n,:));
        %         end
        %         se(g,d) =SE;
        %         lzcomplexity(g,d) = LZC;
        %
        
        
        
        % Uncomment to generate a plot for simulations in all the parameter space
        % ( This is what fills the figure created at the beginning)
        
        %         subplot(length(g),length(d),d+(g-1)*numel(MD))
        %         plot(PSD_G)
        %         xlim ([0 50])
        %         axis off
        
        
        
    end
end

%Save 
save([('Model_Spectral_Features')],'PeakFGlobal', 'PeakFMean','Sync', 'Meta') %,'lzcomplexity','SpecEntropy_Global')

