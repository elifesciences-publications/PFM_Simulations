clear all; close all; clc

clusterthres = 1;
SignalModes = [3 0 20 4 10 25 6 29 26 1 5 16 15 7 23 2 17 41 37 38 47 27 9 8 11 19 21 22 32 24 49 31 14]+1;

% Set paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/')
path(path,'/home/fs0/janineb/scratch/HCP/DMN1200/Functions/')
surfaceL = '/home/fs0/janineb/scratch/HCP/DMN1200/Functions/HCP_S1200_GroupAvg_v1/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii';
surfaceR = '/home/fs0/janineb/scratch/HCP/DMN1200/Functions/HCP_S1200_GroupAvg_v1/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii';
addpath ~steve/NETWORKS/FSLNets;
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions/hline_vline/

% Group maps
in_PFM = '../PFM/maps_full_originalOrder.dtseries.nii'; 
PFMs_group = ft_read_cifti(in_PFM); PFMs_group = PFMs_group.dtseries; PFMs_group(isnan(PFMs_group(:,1))==1,:) = [];

% Create mask
O = ft_read_cifti('../PFM/Results/Spatial_overlap_new.dtseries.nii'); O = O.dtseries(:,1); O(isnan(O)==1,:) = [];
example_cifti = ft_read_cifti(in_PFM);
example_cifti.dtseries = example_cifti.dtseries(:,1); example_cifti.hdr.dim(6) = 1; example_cifti.time = 1;
M = zeros(91282,1); M(O>prctile(O,89.0449)) = 1;
example_cifti.dtseries(isnan(example_cifti.dtseries(:,1))==0,1) = M;
ft_write_cifti('Results/mask',example_cifti,'parameter','dtseries');

% Get cluster
mask = 'Results/mask.dtseries.nii';
PFMs_group = PFMs_group(:,SignalModes);
out_PFMpos = fullfile('Results','ClusterMap_pos.dtseries.nii');
system(sprintf('wb_command -cifti-find-clusters %s %d %d %d %d COLUMN %s -left-surface %s -right-surface %s -cifti-roi %s',in_PFM,clusterthres,100,clusterthres,100,out_PFMpos,surfaceL,surfaceR,mask));
out_PFMneg = fullfile('Results','ClusterMap_neg.dtseries.nii');
system(sprintf('wb_command -cifti-find-clusters %s -%d %d -%d %d COLUMN %s -less-than -left-surface %s -right-surface %s -cifti-roi %s',in_PFM,clusterthres,100,clusterthres,100,out_PFMneg,surfaceL,surfaceR,mask));

% load clustered maps
PFMpos = ft_read_cifti(out_PFMpos); PFMpos = PFMpos.dtseries; 
PFMneg = ft_read_cifti(out_PFMneg); PFMneg = PFMneg.dtseries; 
PFMpos = PFMpos(:,SignalModes); PFMneg = PFMneg(:,SignalModes);
example_cifti = ft_read_cifti(in_PFM);

% LOOP OVER ALL PFMs
masks_group_all = [];
blob_numbers = zeros(size(PFMs_group,2),1);
for n = 1:size(PFMs_group,2)
    
    % extract relevant map
    pos = PFMpos(:,n);
    neg = PFMneg(:,n);
    
    % Create mask with separate clusters
    [~, masks_group] = clusters2rois(pos,neg,[]);
    masks_group(isnan(example_cifti.dtseries(:,1))==1,:) = [];
    masks_group_all = [masks_group_all masks_group];
end
    
masks_group_all = masks_group_all(M==1,:);
masks_group_all(:,sum(masks_group_all)==0) = [];
[C12,munkres_assign] = spatialcorr(masks_group_all,masks_group_all,1);




