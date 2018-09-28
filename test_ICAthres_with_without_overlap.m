clear all; close all; clc

FN = 'LowO_LowM_MedT_NoHighs_01';
dims = [5 7 9 11 13 17 19 21 23 25];

% % Run ICA across different dimensionalities
% for d = dims
%     system(sprintf('fsl_sub melodic -i Results/input_filelist_%s.txt -o Results/Melodic_%s_dim%02d.gica --tr=0.72 --nobet --nomask -a concat --disableMigp -d %d',FN,FN,d, d)); 
% end

% Load data
load(sprintf('Results/PFMsims_atlas_%s.mat',FN),'A','D','atlasParams','params','Pg');
 
% Set paths
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
addpath /vols/Scratch/janineb/HCP/DMN1200/PFM_Simulations/Overlap_functions/

NETcorrs = nan(params.S,15*14/2,3,length(dims));
Tcorrs = nan(params.S,15,3,length(dims));
for d = 1:length(dims)  
    fprintf('dim %d out of %d\n',d, length(dims))
    filename = sprintf('%s_dim%02d',FN,dims(d));
    if dims(d)<15; dims_this = dims(d); else dims_this = 15; end
    params.iN = dims(d);
    [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres_noO,sica_A1_DR_thres_withO] = Melodic_DR(filename,D,atlasParams,params,1);
    
    [C12,munkres_assign] = spatialcorr(sicaPg_new,Pg);
    [i_ica_new,i_Pg] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
    sicaPg_new = sicaPg_new(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1);
    Signs = sign_ica_new; Orders = i_ica_new;
    
    Ai = ones(dims_this); Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);
    
    for s = 1:params.S
        % GT timeseries
        TSG = [A{s}{1}'; A{s}{2}']; TSG = TSG(:,i_Pg);
        NETG = corr(TSG);
        
        % sICA original temporal measures
        a = [sica_A1_DR_new{s}{1}; sica_A1_DR_new{s}{2}];
        r = corr(a(:,Orders).*repmat(Signs',params.T*2,1),a(:,Orders).*repmat(Signs',params.T*2,1));
        r = r-NETG;
        NETcorrs(s,1:An,1,d) = r(Ai);
        r = corr(a(:,Orders).*repmat(Signs',params.T*2,1),TSG); r = r(eye(dims_this)==1);
        Tcorrs(s,1:dims_this,1,d) = r; clear a r
        
        % sICA thresholded temporal measures after removing overlap
        a = [sica_A1_DR_thres_noO{s}{1}; sica_A1_DR_thres_noO{s}{2}];
        r = corr(a(:,Orders).*repmat(Signs',params.T*2,1),a(:,Orders).*repmat(Signs',params.T*2,1));
        r = r-NETG;
        NETcorrs(s,1:An,2,d) = r(Ai);
        r = corr(a(:,Orders).*repmat(Signs',params.T*2,1),TSG); r = r(eye(dims_this)==1);
        Tcorrs(s,1:dims_this,2,d) = r; clear a r
         
        % sICA thresholded temporal measures without removing overlap
        a = [sica_A1_DR_thres_withO{s}{1}; sica_A1_DR_thres_withO{s}{2}];
        r = corr(a(:,Orders).*repmat(Signs',params.T*2,1),a(:,Orders).*repmat(Signs',params.T*2,1));
        r = r-NETG;
        NETcorrs(s,1:An,3,d) = r(Ai);
        r = corr(a(:,Orders).*repmat(Signs',params.T*2,1),TSG); r = r(eye(dims_this)==1);
        Tcorrs(s,1:dims_this,3,d) = r; clear a r
    end
end

T = reshape(NETcorrs,size(NETcorrs,1)*size(NETcorrs,2),3*length(dims));
figure; boxplot(T); set(gcf,'Position',[679 244 999 731],'PaperPositionMode','auto')
vline(3.5:3:30.5,'k')
set(gca,'xtick',2:3:30,'xticklabel',dims)
hline(0,'k')
title('difference between estimated and GT netmats')
xlabel('Dimensionality (GT=15)'); ylabel('edge difference')
print(gcf,'-dpng','-r300','Overlap_tests1.png')

T = reshape(Tcorrs,size(Tcorrs,1)*size(Tcorrs,2),3*length(dims));
figure; boxplot(T); set(gcf,'Position',[679 244 999 731],'PaperPositionMode','auto')
vline(3.5:3:30.5,'k')
title('correlation between estimated and GT timeseries')
set(gca,'xtick',2:3:30,'xticklabel',dims)
xlabel('Dimensionality (GT=15)'); ylabel('Timeseries correlation')
print(gcf,'-dpng','-r300','Overlap_tests2.png')

