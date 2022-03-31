function network_parameter_space(a,MD,expK,C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function loads simulations (Model_Spectral_Features.mat) and create plots
%  1) Input: 
%     
%     a:       bifurcation parameter 
%     MD:      range of delays in ms 
%     expK:    range of couplings (expK)
%     C:       structural connectome, nodes_parcellation_nsubj (i.e. AAL90n32s) 
%
%     Demo:    network_parameter_space(-5,0:1:20,-1:0.1:1.7,'AAL90n32s')
%             
%     Output:  Parameter space plots (Global Peak Frequency, Mean Peak Frequency, Metastability,
%              extra: Spectral Entropy )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Plot GlobalPeakF, MeanPeakF, Meta, Synch

% Uncomment to LOAD Your Simulations: PeakFGlobal, PeakFMean, Sync, Meta parameter range of interest

% if a==-5
%     if C==AAL90n32s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg5')
%         load('Model_Spectral_Features.mat');
%         disp(['Running for' num2str(a) '90AAL32'])
%     elseif C==AAL90n985s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\985AAL\a_neg5')
%         load('Model_Spectral_Features.mat');
%         disp(['Running for' num2str(a) '90AAL985'])
%     elseif C==SHEAF200n32s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[SHEAFER]\Simulations')
%         load('Model_Spectral_Features.mat');
%         disp(['Running for' num2str(a) '200SHEAF32'])
%     end
% 
% elseif a==-0.2
%     if C==AAL90n32s
%     cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg02')
%     load('Model_Spectral_Features.mat');
%     disp(['Running for' num2str(a) '90AAL32'])
%     end
%     
% elseif a==-0.05
%     if C==AAL90n32s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg005\Simulations')
%         load('Model_Spectral_Features.mat');
%         disp(['Running for' num2str(a) '90AAL32'])
%     end
% end

% Load DEMO simulations
if a==-5
    load('Model_Spectral_Features.mat');
    disp(['Running for' num2str(a) '90AAL32'])
end

MD=MD.*1e-3; % Mean Delay in seconds
K=10.^(expK); 

% Create figure
fig=figure ('color', 'w');
colormap(jet);
imagesc(MD*1e3,log10(K),PeakFGlobal)
title('Global Peak Frequency (Hz)','FontSize',20,'FontName','Helvetica')
axis xy
ylabel('Global Coupling K','FontSize',16,'FontName','Helvetica');
xlabel('Mean Delay (ms)','FontSize',16,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
box off
colorbar
caxis([0 45])
h = colorbar;
ylabel(h, 'Frequency (Hz)','FontSize',10,'FontName','Helvetica')
% saveas(fig,'PeakF.png')

% fig=figure ('color', 'w');
% colormap(jet)
% imagesc(MD*1e3,log10(K),PeakFMean)
% title('Mean Peak Frequency (Hz)','FontSize',20,'FontName','Helvetica')
% axis xy
% ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
% xlabel('Mean Delay (ms)','FontSize',14,'FontName','Helvetica')
% yticklabels({'0.1','','1','','10',''});
% box off
% colorbar
% caxis([0 45])
% h = colorbar;
% ylabel(h, 'Frequency (Hz)','FontSize',10,'FontName','Helvetica')
% % saveas(fig,'MeanF.png')

fig=figure ('color', 'w');
colormap(jet)
imagesc(MD*1e3,log10(K),Sync)
axis xy
ylabel('Global Coupling K','FontSize',16,'FontName','Helvetica');
xlabel('Mean Delay (ms)','FontSize',16,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
box off
colorbar
caxis([0 1])
title('Synchrony','FontSize',20,'FontName','Helvetica');
% saveas(fig,'Synch.png')

%
fig=figure ('color', 'w');
colormap(jet)
imagesc(MD*1e3,log10(K),Meta)
title('Metastability','FontSize',20,'FontName','Helvetica');
ylabel('Global Coupling K','FontSize',16,'FontName','Helvetica');
xlabel('Mean Delay (ms)','FontSize',16,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
box off
axis xy
colorbar
% caxis([0 0.2])
% saveas(fig,'Meta.png')
% subplot(2,3,3)
% imagesc(K,MD*1e3,PeakF')
% axis xy
% ylabel('Mean Delay (ms)')
% xlabel('Coupling K')
% colorbar
% title('Oscillation Frequency')
% colormap(jet)

%% If you also measured Spectral Entropy in Function 1 
%  in DEMO simulations there is no Spectral Entropy calculation

% fig=figure ('color', 'w');
% colormap(jet)
% imagesc(MD*1e3,log10(K),SpecEntropy_Global)
% title('Spectral Entropy','FontSize',20,'FontName','Helvetica');
% ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica');
% xlabel('Mean Delay (ms)','FontSize',14,'FontName','Helvetica')
% yticklabels({'0.1','','1','','10',''});
% axis xy
% colorbar
% caxis([0 1])
% % saveas(fig,'Spec_Entropy.png')

% subplot(2,3,3)
% imagesc(K,MD*1e3,PeakF')
% axis xy
% ylabel('Mean Delay (ms)')
% xlabel('Coupling K')
% colorbar
% title('Oscillation Frequency')
% colormap(jet)


