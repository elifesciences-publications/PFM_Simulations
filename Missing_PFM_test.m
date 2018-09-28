clear all; close all; clc

Ins = {'LowO_LowM_MedT_NoHighs'};
warning off

% Set paths
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
addpath /vols/Scratch/janineb/HCP/DMN1200/PFM_Simulations/Overlap_functions/

% Load data
load(sprintf('Results/PFMsims_atlas_%s_01.mat',Ins{1}),'atlasParams','params');

MISSING = zeros(100,6);
I = 1;

for d = 1:params.nRepeats
     
    filename = sprintf('%s_%02d',Ins{1},d);
    
    % Load data
    load(sprintf('Results/PFMsims_atlas_%s.mat',filename),...
        'Pg','P','A','D','atlasParams','params');
    fprintf('Running iteration %d out of %d\n',d,params.nRepeats);
    
    % Rerun PFM gathering
    [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(filename,params);
    
    % Get all group maps and fix order & sign to match ground truth:
    pfmPg_new = pfmPg1_new;
    [C12,munkres_assign] = spatialcorr(pfmPg1_new,Pg);
    [i_pfm_new,~] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
    pfmPg_new = pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
    
    % Get group nets
    NET = zeros(params.iN,params.iN,params.S);
    maps = zeros(atlasParams.V,params.S,params.iN);
    maps_corrs = zeros(params.S*(params.S-1)/2,params.iN);
    Ai = ones(params.S); Ai = triu(Ai,1); Ai = find(Ai==1);
    for s = 1:params.S 
        R = corr([A{s}{1}'; A{s}{2}']); R = triu(R,1);
        NET(:,:,s) = R;  
        M = P{s};
        M = M(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
        maps(:,s,:) = M;
    end
    for n = 1:params.iN
        r = corr(squeeze(maps(:,:,n)));
        maps_corrs(:,n) = r(Ai);
    end
    maps_corrs = mean(maps_corrs);
    
    % Find missing group maps
    [C12,munkres_assign] = spatialcorr(pfmPg_new,Pg);
    C = C12(eye(params.iN)==1);
    MIS = find(C<0.5);
    MISSING(I:I+length(MIS)-1,1) = C(MIS);
    for n = 1:length(MIS)
        [~,i] = find(max(abs(C12(:,MIS(n)))));
        MISSING(I+n-1,2) = C12(i,MIS(n));
        MISSING(I+n-1,3) = maps_corrs(MIS(n));
        MISSING(I+n-1,4) = mean(maps_corrs);
        MISSING(I+n-1,5) = max(squeeze(abs(NET(MIS(n),i,:))));
        MISSING(I+n-1,6) = max(abs(NET(:)));
    end
    I = I+length(MIS);
    
end

