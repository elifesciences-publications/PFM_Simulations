clear all; close all; clc

Ins = {'original','HighO_LowM_LowT','LowO_LowM_HighT','LowO_HighM_LowT'};
n=3;

% Load stuff
load(sprintf('Results/PFMsims_atlas_%s.mat',Ins{n}),...
'pfmPg1_new','Pg','ticaP1','sicaPg_new','P','A',...
'pfmP1_new','pfmA1_new','ticaP1_DR','ticaA1_DR',...
'sica_P1_DR_new','sica_A1_DR_new','D','atlasParams','params','P');

% Run stuff
[sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres] = Melodic_DR(Ins{n},D,atlasParams,params);
[pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(Ins{n},params);

pfmPg_new = pfmPg1_new;
[C12,munkres_assign] = spatialcorr(pfmPg_new,Pg);
[i_pfm_new,~] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
pfmPg_new = pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
GroupMaps(:,:,2,n) = pfmPg_new; Signs(:,2,n) = sign_pfm_new; Orders(:,2,n) = i_pfm_new;
ticaPg = ticaP1{1};
[C12,munkres_assign] = spatialcorr(ticaPg,Pg);
[i_tica,~] = find(munkres_assign==1); sign_tica = sign(C12(munkres_assign==1));
ticaPg = ticaPg(:,i_tica).*repmat(sign_tica',atlasParams.V,1);
GroupMaps(:,:,3,n) = ticaPg; Signs(:,3,n) = sign_tica; Orders(:,3,n) = i_tica;
[C12,munkres_assign] = spatialcorr(sicaPg_new,Pg);
[i_ica_new,~] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
sicaPg_new = sicaPg_new(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1);
GroupMaps(:,:,4,n) = ticaPg; Signs(:,4,n) = sign_ica_new; Orders(:,4,n) = i_ica_new;

%%%%% Comparison time:
for s = 1
    
    testP1_ica = sica_P1_DR_new{s}(:,Orders(:,4,n)).*repmat(sign_ica_new',atlasParams.V,1);
    testP1_pfm = pfmP1_new{s}(:,Orders(:,2,n)).*repmat(sign_pfm_new',atlasParams.V,1);
    testP1_gt = P{s};
    
    % Set max to one for better visualisation
    testP1_ica_norm = testP1_ica.*repmat(sign(mean(testP1_gt)),size(testP1_ica,1),1)./repmat(max(abs(testP1_ica)),size(testP1_ica,1),1);
    testP1_pfm_norm = testP1_pfm.*repmat(sign(mean(testP1_gt)),size(testP1_pfm,1),1)./repmat(max(abs(testP1_pfm)),size(testP1_pfm,1),1);
    testP1_gt_norm = testP1_gt.*repmat(sign(mean(testP1_gt)),size(testP1_gt,1),1)./repmat(max(abs(testP1_gt)),size(testP1_gt,1),1);
    
    % spatial nets:
    %figure; imagesc(corr(testP1_pfm),[-1 1])
    %figure; imagesc(corr(testP1_ica),[-1 1])
    %figure; imagesc(corr(testP1_gt),[-1 1])
    
    % Pick an interesting edge (2-13), which is positive in GT & PFM, but negative in ICA
    figure; set(gcf,'Position',[10 10 750 955],'PaperPositionMode','auto')
    imagesc([testP1_gt_norm(:,[2 13]) testP1_pfm_norm(:,[2 13]) testP1_ica_norm(:,[2 13])])
    vline(2.5:2:6.5,'k'); title('Ground Truth                              PFM                                   ICA');
    set(gca,'xtick',1:6,'xticklabel',{'GT2','GT13','PFM2','PFM13','ICA2','ICA13'})
    mapcorrs = corr([testP1_gt(:,[2 13]) testP1_pfm(:,[2 13]) testP1_ica(:,[2 13])])
end
print(gcf,'-dpng','-r300','Results/Comparison_example_subject.png')


% Set max to one for better visualisation
testPg_ica_norm = sicaPg_new.*repmat(sign(mean(Pg)),size(sicaPg_new,1),1)./repmat(max(abs(sicaPg_new)),size(sicaPg_new,1),1);
testPg_pfm_norm = pfmPg_new.*repmat(sign(mean(Pg)),size(pfmPg_new,1),1)./repmat(max(abs(pfmPg_new)),size(pfmPg_new,1),1);
testPg_gt_norm = Pg.*repmat(sign(mean(Pg)),size(Pg,1),1)./repmat(max(abs(Pg)),size(Pg,1),1);

figure; set(gcf,'Position',[760 10 750 955],'PaperPositionMode','auto')
imagesc([testPg_gt_norm(:,[2 13]) testPg_pfm_norm(:,[2 13]) testPg_ica_norm(:,[2 13])])
vline(2.5:2:6.5,'k'); title('Ground Truth                              PFM                                   ICA');
set(gca,'xtick',1:6,'xticklabel',{'GT2','GT13','PFM2','PFM13','ICA2','ICA13'})
print(gcf,'-dpng','-r300','Results/Comparison_example_group.png')




