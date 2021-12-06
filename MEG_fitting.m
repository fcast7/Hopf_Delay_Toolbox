function MEG_fitting 

%load DEMO
addpath('Hopf_Delay_Toolbox');
load('MEG_Fitting.mat'); 
expK=-1:0.1:1.7;
K=10.^(expK);

%%

for s=1:89
    Dist_MEG= squeeze(Error_MEG_PSD_Sub(:,:,s));
%     Dist_MEG= squeeze(Error_MEG_PSD_Sub(:,:,s));
    [~, index_best_fit]= min(Dist_MEG(:));
    [k_s(s),D_s(s)]=ind2sub(size(Dist_MEG),index_best_fit);
%     D_s(s)=find(sum(squeeze(Corr_MEG_PSD_Sub(:,:,s))==max(max(squeeze(Corr_MEG_PSD_Sub(:,:,s))))));
%     k_s(s)=find(sum(squeeze(Corr_MEG_PSD_Sub(:,:,s))==max(max(squeeze(Corr_MEG_PSD_Sub(:,:,s)))),2));
    disp(['Best fit to sub' num2str(s) ' for delay= ' num2str(MD(D_s(s))*1000) 'ms and k=1E' num2str(expK(k_s(s)))])
end

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
plot(MD(D_s)*1000,expK(k_s),'*w');
legend('Best fit of individual MEG PSD')


for s=1:89
    Corr_MEG= squeeze(Corr_MEG_PSD_Sub(:,:,s));
    [~, index_best_fit]= max(Corr_MEG(:));
    [k_s(s),D_s(s)]=ind2sub(size(Corr_MEG),index_best_fit);
%     D_s(s)=find(sum(squeeze(Corr_MEG_PSD_Sub(:,:,s))==max(max(squeeze(Corr_MEG_PSD_Sub(:,:,s))))));
%     k_s(s)=find(sum(squeeze(Corr_MEG_PSD_Sub(:,:,s))==max(max(squeeze(Corr_MEG_PSD_Sub(:,:,s)))),2));
    disp(['Best fit to sub' num2str(s) ' for delay= ' num2str(MD(D_s(s))*1000) 'ms and k=1E' num2str(expK(k_s(s)))])
end

figure ('color','w')
colormap(jet)
imagesc(MD*(1E3),log10(K),Fit_MEG_PSD) 
title('Fit MEG PSD','FontSize',14,'FontName','Helvetica')
ylabel('Global Coupling K','FontSize',14,'FontName','Helvetica')
xlabel('Mean delay (ms)','FontSize',14,'FontName','Helvetica')
yticklabels({'0.1','','1','','10',''});
axis xy
colorbar
h = colorbar;
ylabel(h, 'Correlation','FontSize',11,'FontName','Helvetica')
hold on
plot(MD(D_s)*1000,expK(k_s),'*w');
legend('Best fit of individual MEG PSD')




