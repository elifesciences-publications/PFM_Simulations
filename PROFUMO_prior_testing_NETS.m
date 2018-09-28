clear all; close all; clc

% Rerun PFM gathering
Ins = {'MedO_MedM_MedT_01_strong20','MedO_MedM_MedT_01_med20','MedO_MedM_MedT_01_weak20'};

% Load data
load(sprintf('Results/PFMsims_atlas_%s.mat','MedO_MedM_MedT_01'),...
    'Pg','P','A','D','atlasParams','params');

% GT group maps
GroupMaps(:,:,length(Ins)+1) = Pg;

% Initiate figures
Sfig = figure; set(Sfig,'Position',[0 0 1950 950],'PaperPositionMode','auto')
Tfig = figure; set(Tfig,'Position',[0 0 1950 950],'PaperPositionMode','auto')

% GT netmats
for s = 1:10%params.S
    [r,~] = spatialcorr(P{s},P{s}); 
    set(0,'CurrentFigure',Sfig); subplot(length(Ins)+1,10,s); imagesc(r,[-0.5 0.5]);
    title(sprintf('Subject %d Snets',s)); ylabel('GT')
    r = corr([A{s}{1}'; A{s}{2}']);
    set(0,'CurrentFigure',Tfig); subplot(length(Ins)+1,10,s); imagesc(r,[-0.5 0.5]);
    title(sprintf('Subject %d Tnets',s)); ylabel('GT')
end

for n = 1:length(Ins)
    [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(Ins{n},params);
    
    pfmPg_new = pfmPg1_new;
    [C12,munkres_assign] = spatialcorr(pfmPg1_new,Pg);
    [i_pfm_new,~] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
    pfmPg_new = pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
    GroupMaps(:,:,n) = pfmPg_new; Signs(:,2,n) = sign_pfm_new; Orders(:,2,n) = i_pfm_new;
    
    for s = 1:10%params.S
        
        % PFM nets
        p = pfmP1_new{s}; [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1), p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1));
        set(0,'CurrentFigure',Sfig); subplot(length(Ins)+1,10,s+(n*10)); imagesc(r,[-0.5 0.5]);
        if s==1; ylabel(sprintf('%s',Ins{s}),'interpreter','none'); end
        [~,a] = cov2corr(inv(pfmNet{s}));
        a = a(:,Orders(:,2,n)); a = a(Orders(:,2,n),:);
        a = a.*(repmat(Signs(:,2,n)',params.N,1) .* repmat(Signs(:,2,n),1,params.N));
        set(0,'CurrentFigure',Tfig); subplot(length(Ins)+1,10,s+(n*10)); imagesc(a,[-0.5 0.5]);
       if s==1; ylabel(sprintf('%s',Ins{s}),'interpreter','none'); end
    end
end

