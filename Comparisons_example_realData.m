clear all; close all; clc
% Indices (matched maps based on Spatial_correlation_ICA.m):
Iica = [1 23 2 20 4 3 5 8 11 19 10];
Ipfm = [1 2 12 3 7 11 22 10 25 15 18];

addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
HCPdir = '/vols/Data/HCP/Phase2/subjects1200';
INdir = '/vols/Scratch/HCP/rfMRI/PROFUMO/S1003_M50.pfm/FinalModel';
scans = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
Pdir = '/vols/Scratch/HCP/rfMRI/PROFUMO/S1003_M50.pfm/FinalModel/';
subs = dir('/vols/Scratch/HCP/rfMRI/PROFUMO/S820_M50_Aug16.pfm/FinalModel/Subjects/');
subs = subs(3:end); subs(169) = [];

% Group maps:
groupICA = '/vols/Data/HCP/Phase2/group1200/groupICA/groupICA_3T_HCP1200_MSMAll_d25.ica/melodic_IC.dtseries.nii';
gICA = ft_read_cifti(groupICA); gICA = gICA.dtseries(:,Iica); gICA(isnan(gICA(:,1))==1,:) = [];
groupPFM = '../PFM/maps_full_originalOrder_max.dtseries.nii';
gPFM = ft_read_cifti(groupPFM); gPFM = gPFM.dtseries(:,Ipfm); gPFM(isnan(gPFM(:,1))==1,:) = [];
gICA = gICA.*repmat(sign(mean(gPFM)),size(gICA,1),1)./repmat(max(abs(gICA)),size(gICA,1),1);
gPFM = gPFM.*repmat(sign(mean(gPFM)),size(gPFM,1),1)./repmat(max(abs(gPFM)),size(gPFM,1),1);

% Edge to look at:
Ns = [3 7]; %[5 1] [10 4] [3 7]; 3==DMN | 4==POS2

% Sort into suitable order
[Ms1,Is1] = sort(gPFM(:,Ns(2)),'descend');
[Ms2,Is2] = sort(gPFM(:,Ns(1)),'descend');
Is = [Is1(1:find(Ms1<0.2,1,'first'))' Is2(1:find(Ms2<0.2,1,'first'))']; Is = [Is setdiff(1:91282,Is)];
gICA = gICA(Is,:); gPFM = gPFM(Is,:);
Grays = 1:find(Ms1<0.2,1,'first')+find(Ms2<0.2,1,'first');

for s=6;
    means = h5read(fullfile(Pdir,'Subjects',subs(s).name,'SpatialMaps.post','Signal','Means.hdf5'),'/dataset');
    probs = h5read(fullfile(Pdir,'Subjects',subs(s).name,'SpatialMaps.post','Signal','MembershipProbabilities.hdf5'),'/dataset');
    maps = means .* probs; maps = maps(:,Ipfm); clear means probs
    
    icamaps = ft_read_cifti(sprintf('/vols/Data/HCP/Phase2/group1200/node_maps/3T_HCP1200_MSMAll_d25_ts2_Z/%s.dtseries.nii',subs(s).name));
    icamaps = icamaps.dtseries; icamaps(isnan(icamaps(:,1))==1,:) = []; icamaps = icamaps(:,Iica);
    
    maps_norm = maps.*repmat(sign(mean(maps)),size(maps,1),1)./repmat(max(abs(maps)),size(maps,1),1);
    icamaps_norm = icamaps.*repmat(sign(mean(maps)),size(icamaps,1),1)./repmat(max(abs(icamaps)),size(icamaps,1),1);
    
    icamaps_norm = icamaps_norm(Is,:); maps_norm = maps_norm(Is,:);
    
    figure; set(gcf,'Position',[0 10 750 955],'PaperPositionMode','auto')
    imagesc([maps_norm(Grays,Ns) icamaps_norm(Grays,Ns) gPFM(Grays,Ns) gICA(Grays,Ns)],[-1 1])
    vline(2.5:2:8.5,'k');
    Names = {sprintf('PFM%d',Ns(1)),sprintf('PFM%d',Ns(2)),sprintf('ICA%d',Ns(1)),sprintf('ICA%d',Ns(2)),...
        sprintf('groupPFM%d',Ns(1)),sprintf('groupPFM%d',Ns(2)),sprintf('groupICA%d',Ns(1)),sprintf('groupICA%d',Ns(2))};
    set(gca,'xtick',1:8,'xticklabel',Names);
    print(gcf,'-dpng','-r300','Results/Comparison_example_data.png')
    
    maps_thres = maps_norm; maps_thres(abs(maps_norm)<0.2) = 0;
    icamaps_thres = icamaps_norm; icamaps_thres(abs(icamaps_norm)<0.2) = 0;
    gPFM_thres = gPFM; gPFM_thres(abs(gPFM)<0.2) = 0;
    gICA_thres = gICA; gICA_thres(abs(gICA)<0.2) = 0;
    
%     figure; 
%     subplot(2,1,1); 
%     imagesc(corr([maps_norm(:,Ns) icamaps_norm(:,Ns) gPFM(:,Ns) gICA(:,Ns)]),[-1 1])
%     hline(2.5:2:8.5,'k'); vline(2.5:2:8.5,'k');
%     title('spatial correlations no thresholding'); set(gca,'xtick',1:8,'xticklabel',Names,'ytick',1:8,'yticklabel',Names);
%     subplot(2,1,2); 
%     imagesc(corr([maps_thres(:,Ns) icamaps_thres(:,Ns) gPFM(:,Ns) gICA(:,Ns)]),[-1 1])
%     hline(2.5:2:8.5,'k'); vline(2.5:2:8.5,'k');
%     title('spatial correlations no thresholding'); set(gca,'xtick',1:8,'xticklabel',Names,'ytick',1:8,'yticklabel',Names);
    
end