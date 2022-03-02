function meg_psd_model_screening(a, MD, expK, C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function loads simulations (MEG_sensor_PSD_Fitting.mat) and geenrates plots
%  1) Input: 
%     
%     a:       bifurcation parameter 
%     MD:      range of delays in ms 
%     expK:    range of couplings (expK)
%     C:       structural connectome, nodes_parcellation_nsubj (i.e. AAL90n32s) 
%
%     Example: network_parameter_space(-0.05,0:1:20,-1:0.1:1.7,AAL90n32s)
%     Demo:    network_parameter_space(-5,0:1:20,-1:0.1:1.7,AAL90n32s)
%             
%     Output:  Model disparity plot (Squared Euclidean Distance). Best fitting 
%              points visualised both for each subject and their mean
%
% Example: meg_psd_model_screening(-5,0:1:20,-1:0.1:1.7, AAL90n32s)
% Code by Francesca Castaldo, 2022 francesca.castaldo.20@ucl.ac.uk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add the toolbox to your own path
addpath(genpath('C:\Users\fcast\OneDrive - University College London\Cabral_Castaldo\Code\Hopf_Delay_Toolbox\Hopf_Delay_Toolbox_2ndVersion'))

% if a==-5
%     if C==AAL90n32s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg5')
%         load ('MEG_sensor_PSD_Fitting.mat');
%         K=10.^(expK);
%         disp(['Running for' num2str(a) '90AAL32'])
%     elseif C==AAL90n985s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\985AAL\a_neg5')
%         load ('MEG_sensor_PSD_Fitting.mat');
%         K=10.^(expK);
%         disp(['Running for' num2str(a) '90AAL985'])
%     elseif C==SHEAF200n32s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[SHEAFER]\Simulations')
%         load ('MEG_sensor_PSD_Fitting.mat');
%         K=10.^(expK);
%         disp(['Running for' num2str(a) '200SHEAF32'])
%     end
% 
% elseif a==-0.2
%     if C==AAL90n32s
%     cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg02')
%     load ('MEG_sensor_PSD_Fitting.mat');
%     K=10.^(expK);
%     disp(['Running for' num2str(a) '90AAL32'])
%     end
%     
% elseif a==-0.05
%     if C==AAL90n32s
%         cd('C:\Users\fcast\OneDrive - University College London\CLUSTER\PROJECT[AAL]\New_Simu\32AAL\a_neg005\Simulations')
%         load ('MEG_sensor_PSD_Fitting.mat');
%         K=10.^(expK);
%         disp(['Running for' num2str(a) '90AAL32'])
%     end
% end


% Demo Data:
if a==-5 && C==AAL90n32s
    load ('MEG_sensor_PSD_Fitting.mat');
end


%% Fit for each subject

for s=1:89
    Dist_MEG= squeeze(Error_MEG_PSD_Sub(:,:,s));
    [~, index_best_fit]= min(Dist_MEG(:));
    best_fit_subj(s)=index_best_fit;
    [k_s(s),D_s(s)]=ind2sub(size(Dist_MEG),index_best_fit);
    disp(['Best fit to sub' num2str(s) ' for delay= ' num2str(MD(D_s(s))*1000) 'ms and k=1E' num2str(expK(k_s(s))) ' index ' num2str(index_best_fit)])
end

x=unique(best_fit_subj);
N=numel(x);

count = zeros(N,1);
for ii = 1:N
    count(ii) = sum(best_fit_subj==x(ii));
end

disp([ x(:) count ]);


figure ('color','w')
colormap(jet)
imagesc(MD*1E3,log10(K),Error_MEG_PSD)%,'AlphaData',~isnan(Error_MEG_PSD)
title('Distance MEG PSD','FontSize',14,'FontName','Helvetica')
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica')
xlabel('Mean delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
axis xy
colorbar
h = colorbar;
ylabel(h, 'Squared Euclidean Distance','FontSize',11,'FontName','Helvetica')


hold on
% plot(MD(D_s)*1000,expK(k_s),'*w');
for s=1:89
plot(MD(D_s(s))*1000,expK(k_s(s)),'*w');
legend('Best fit of individual MEG PSD')
end


%% Fit the subject mean 

Dist_MEG_mean=Error_MEG_PSD;
[~, index_best_fit_mean]= min(Dist_MEG_mean(:));
%     best_fit_meansubj(g,d)=index_best_fit_mean;
[k_s_mean,D_s_mean]=ind2sub(size(Dist_MEG_mean),index_best_fit_mean);
disp(['Best fit for delay= ' num2str(MD(D_s_mean)*1000) 'ms and k=1E' num2str(expK(k_s_mean)) ' index ' num2str(index_best_fit_mean)])

figure ('color','w')
colormap(jet)
imagesc(MD*1E3,log10(K),Error_MEG_PSD)%,'AlphaData',~isnan(Error_MEG_PSD)
title('Distance MEG PSD','FontSize',14,'FontName','Helvetica')
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica')
xlabel('Mean delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
axis xy
colorbar
h = colorbar;
ylabel(h, 'Squared Euclidean Distance','FontSize',11,'FontName','Helvetica')

hold on
plot(MD(D_s_mean)*1000,expK(k_s_mean),'*w');
legend('Best fit of mean MEG PSD')


%% Correlation

% for s=1:89
%     Corr_MEG= squeeze(Corr_MEG_PSD_Sub(:,:,s));
%     [~, index_best_fit]= max(Corr_MEG(:));
%     
%     best_fit_subj(s)=index_best_fit;
%     [k_s(s),D_s(s)]=ind2sub(size(Corr_MEG),index_best_fit);
% %     D_s(s)=find(sum(squeeze(Corr_MEG_PSD_Sub(:,:,s))==max(max(squeeze(Corr_MEG_PSD_Sub(:,:,s))))));
% %     k_s(s)=find(sum(squeeze(Corr_MEG_PSD_Sub(:,:,s))==max(max(squeeze(Corr_MEG_PSD_Sub(:,:,s)))),2));
%     disp(['Best fit to sub' num2str(s) ' for delay= ' num2str(MD(D_s(s))*1000) 'ms and k=1E' num2str(expK(k_s(s)))])
% end

% x=unique(best_fit_subj);
% N=numel(x);
% 
% count = zeros(N,1);
% for ii = 1:N
%     count(ii) = sum(best_fit_subj==x(ii));
% end
% 
% disp([ x(:) count ]);




% figure ('color','w')
% colormap(jet)
% imagesc(MD*(1E3),log10(K),Fit_MEG_PSD) 
% title('Fit MEG PSD','FontSize',14,'FontName','Helvetica')
% ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica')
% xlabel('Mean delay (ms)','FontSize',14,'FontName','Helvetica')
% yticklabels({'0.1','','1','','10',''});
% axis xy
% colorbar
% h = colorbar;
% ylabel(h, 'Correlation','FontSize',11,'FontName','Helvetica')
% hold on
% 
% plot(MD(D_s)*1000,expK(k_s),'*w');
% legend('Best fit of individual MEG PSD', 'Location', 'northoutside')
% 

