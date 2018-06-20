clear all; close all; clc

% Only works if simulations have the same number of voxels, subjects, and modes

% Inputs to compare
Ins = {'original','HighO_LowM_LowT','LowO_LowM_HighT','LowO_HighM_LowT'};
PFMnet = 'TS'; % 'PROFUMO' or 'TS' or 'DR'
ICA_Snet = 'original'; % 'thres_Smap' or 'original'
ICA_Tnet = 'original'; % 'thres_Gmap' or 'original'
tICA = 'high_dimICA'; % 'high_dimICA' or 'low_dimPCA'
    
% Set paths
warning off
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
addpath /vols/Scratch/janineb/HCP/DMN1200/PFM_Simulations/Overlap_functions/

% Intitialize variables
load(sprintf('Results/PFMsims_atlas_%s.mat',Ins{1}),'atlasParams','params');
Ai = ones(params.iN); Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);
GroupMaps = zeros(atlasParams.V,params.iN,4,length(Ins));
Signs = zeros(params.iN,4,length(Ins));
Orders = zeros(params.iN,4,length(Ins));
Snets = zeros(params.S,An,4,length(Ins));
Tnets = zeros(params.S,An,4,length(Ins));
Scorrs = zeros(params.S,params.iN,4,length(Ins));
AllNames = {'Ground truth','PFM','tICA','sICA'};
        
%%% Loop over all comparison datasets and extract relevant information:
for n = 1:length(Ins)
    fprintf('Loading data for input %d (%s)\n',n,Ins{n})
    
    % Load data
    load(sprintf('Results/PFMsims_atlas_%s.mat',Ins{n}),...
        'pfmPg1_new','Pg','ticaP1','P','A','sticaP1_DR','sticaA1_DR',...
        'pfmP1_new','pfmA1_new','ticaP1_DR','ticaA1_DR',...
        'D','atlasParams','params','P'); 
    
    if strcmp(tICA,'high_dimICA');
        ticaP1_DR = sticaP1_DR;
        ticaA1_DR = sticaA1_DR;
    end
    
    % Rerun DR
    [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres] = Melodic_DR(Ins{n},D,atlasParams,params);
    
    % Rerun PFM gathering
    [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(Ins{n},params); 
    
    % Get all group maps and fix order & sign to match ground truth:
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
    
    GroupMaps(:,:,1,n) = Pg;
    clear pfmPg_new C12 munkres_assign i_pfm_new sign_pfm_new ticaPg i_tica sign_tica i_ica_new sign_ica_new
    
    % Get subject temporal and spatial netmats
    Nedge = params.iN*params.iN;
    for s = 1:params.S
        
        % Ground truth nets
        [r,~] = spatialcorr(P{s},P{s}); Snets(s,:,1,n) = r(Ai);
        r = corr([A{s}{1}'; A{s}{2}']); Tnets(s,:,1,n) = r(Ai);
        
        % PFM nets
        p = pfmP1_new{s}; [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1), p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1)); 
        Snets(s,:,2,n) = r(Ai);
        [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1),P{s}); 
        Scorrs(s,:,2,n) = r(eye(params.iN)==1)'; 
        if strcmp(PFMnet,'TS')
            a = [pfmA1_new{s}{1}'; pfmA1_new{s}{2}']; r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1)); 
            Tnets(s,:,2,n) = r(Ai); 
        elseif strcmp(PFMnet,'PROFUMO')
            [~,a] = cov2corr(inv(pfmNet{s})); 
            a = a(:,Orders(:,2,n)); a = a(Orders(:,2,n),:);
            a = a.*(repmat(Signs(:,2,n)',params.iN,1) .* repmat(Signs(:,2,n),1,params.iN));
            Tnets(s,:,2,n) = a(Ai); 
        elseif strcmp(PFMnet,'DR')
            a = [];
            for r = 1:params.R(1)
                a = [a; nets_demean((pinv(nets_demean(p)) * nets_demean(double(D{s}{r})))')];
            end
            r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1)); 
            Tnets(s,:,2,n) = r(Ai); 
        end
        clear p a
        
        % tICA nets
        p = ticaP1_DR{s}; [r,~] = spatialcorr(p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1), p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1)); 
        Snets(s,:,3,n) = r(Ai); 
        [r,~] = spatialcorr(p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1),P{s}); 
        Scorrs(s,:,3,n) = r(eye(params.iN)==1); clear p
        a = [ticaA1_DR{s}{1}'; ticaA1_DR{s}{2}']; r = corr(a(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',params.T*2,1)); 
        Tnets(s,:,3,n) = r(Ai); clear a
        
        % sICA nets
        p = sica_P1_DR_thres{s}; 
        if strcmp(ICA_Snet,'original')
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1), p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1)); 
            Snets(s,:,4,n) = r(Ai); 
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1),P{s}); 
            Scorrs(s,:,4,n) = r(eye(params.iN)==1);
        elseif strcmp(ICA_Snet,'thres_Smap')
            p(abs(p)<5) = 0; 
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1), p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1)); 
            Snets(s,:,4,n) = r(Ai); 
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1),P{s}); 
            Scorrs(s,:,4,n) = r(eye(params.iN)==1); clear p
        end
        if strcmp(ICA_Tnet,'original')          
            a = [sica_A1_DR_new{s}{1}; sica_A1_DR_new{s}{2}]; r = corr(a(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',params.T*2,1)); 
            Tnets(s,:,4,n) = r(Ai); clear a
        elseif strcmp(ICA_Tnet,'thres_Gmap')          
            a = [sica_A1_DR_thres{s}{1}; sica_A1_DR_thres{s}{2}]; r = corr(a(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',params.T*2,1)); 
            Tnets(s,:,4,n) = r(Ai); clear a
        end
    end
    clear('pfmPg1_new','Pg','ticaP1','sicaPg_new','P','A',...
        'pfmP1_new','pfmA1_new','ticaP1_DR','ticaA1_DR',...
        'sica_P1_DR_new','sica_A1_DR_new','D','P'); 
end

% Plot comparison of group maps
figure; set(gcf,'Position',[0 570 1900 955],'PaperPositionMode','auto')
for n = 1:length(Ins)
    subplot(length(Ins),4,(n-1)*4+1); imagesc(GroupMaps(:,:,1,n)); colormap parula; set(gca,'xtick',1:params.iN); grid on
    if n == 1; title(sprintf('%s group maps',AllNames{1}),'color','k'); end
    ylabel(sprintf('%s',Ins{n}),'interpreter','none')
    
    subplot(length(Ins),4,(n-1)*4+2); imagesc(GroupMaps(:,:,2,n)); colormap parula; set(gca,'xtick',1:params.iN); grid on
    for i = 1:params.iN; r = corr(GroupMaps(:,i,1,n),GroupMaps(:,i,2,n)); h = text(i-0.1,atlasParams.V-0.05*atlasParams.V,sprintf('%1.2f',r)); set(h,'rotation',90); end
    if n == 1; title(sprintf('%s group maps',AllNames{2}),'color','b'); end
    
    subplot(length(Ins),4,(n-1)*4+3); imagesc(GroupMaps(:,:,3,n)); colormap parula; set(gca,'xtick',1:params.iN); grid on
    for i = 1:params.iN; r = corr(GroupMaps(:,i,1,n),GroupMaps(:,i,3,n)); h = text(i-0.1,atlasParams.V-0.05*atlasParams.V,sprintf('%1.2f',r)); set(h,'rotation',90); end
    if n == 1; title(sprintf('%s group maps',AllNames{3}),'color','r'); end
    
    subplot(length(Ins),4,(n-1)*4+4); imagesc(GroupMaps(:,:,4,n)); colormap parula; set(gca,'xtick',1:params.iN); grid on
    for i = 1:params.iN; r = corr(GroupMaps(:,i,1,n),GroupMaps(:,i,4,n)); h = text(i-0.1,atlasParams.V-0.05*atlasParams.V,sprintf('%1.2f',r)); set(h,'rotation',90); end
    if n == 1; title(sprintf('%s group maps',AllNames{4}),'color',[0.6 0 1]); end
end
print(gcf,'-dpng','-r300','Results/Comparison_GroupMaps.png')

% Plot comparison of subject simlarity of temporal netmats
figure; set(gcf,'Position',[0 0 1550 955],'PaperPositionMode','auto')
for n = 1:length(Ins)
    subplot(length(Ins),3,(n-1)*3+1);
    r = corr(Tnets(:,:,1,n)',Tnets(:,:,2,n)'); imagesc(r,[-1 1]); colorbar; colormap parula
    if n == 1; title(sprintf('Tnets %s (%s)',AllNames{2},PFMnet),'color','b','interpreter','none'); end
    r = r(eye(params.S)==1); xlabel(sprintf('mean r=%1.2f',mean(r)));
    ylabel(sprintf('%s',Ins{n}),'interpreter','none')
    
    subplot(length(Ins),3,(n-1)*3+2);
    r = corr(Tnets(:,:,1,n)',Tnets(:,:,3,n)'); imagesc(r,[-1 1]); colorbar; colormap parula
    if n == 1; title(sprintf('Tnets %s',AllNames{3}),'color',[0.6 0 1],'interpreter','none'); end
    r = r(eye(params.S)==1); xlabel(sprintf('mean r=%1.2f',mean(r)));
    
    subplot(length(Ins),3,(n-1)*3+3);
    r = corr(Tnets(:,:,1,n)',Tnets(:,:,4,n)'); imagesc(r,[-1 1]); colorbar; colormap parula
    if n == 1; title(sprintf('Tnets %s (%s)',AllNames{4},ICA_Tnet),'color','r','interpreter','none'); end
    r = r(eye(params.S)==1); xlabel(sprintf('mean r=%1.2f',mean(r)));   
end
print(gcf,'-dpng','-r300','Results/Comparison_TemporalNets.png')
  
% Plot temporal netmat scatter plots 
figure; set(gcf,'Position',[0 0 1900 955],'PaperPositionMode','auto')
for n = 1:length(Ins)
    Cg = reshape(Tnets(:,:,1,n),params.S*An,1);
    
    subplot(length(Ins),3,(n-1)*3+1);
    Cn = reshape(Tnets(:,:,2,n),params.S*An,1); r = corr(Cg,Cn);
    scatplot(Cg,Cn); axis([-1 1 -1 1]); hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else plot(-1:0.1:1,1:-0.1:-1,'r'); end
    if n == 1; title(sprintf('Tnets scatterplot: %s (%s)',AllNames{2},PFMnet),'color','b','interpreter','none'); end
    xlabel(sprintf('Ground Truth: r=%1.2f',r));
    ylabel(sprintf('%s',Ins{n}),'interpreter','none')
    
    subplot(length(Ins),3,(n-1)*3+2);
    Cn = reshape(Tnets(:,:,3,n),params.S*An,1); r = corr(Cg,Cn);
    scatplot(Cg,Cn); axis([-1 1 -1 1]); hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else plot(-1:0.1:1,1:-0.1:-1,'r'); end
    if n == 1; title(sprintf('Tnets scatterplot: %s',AllNames{3}),'color',[0.6 0 1],'interpreter','none'); end
    xlabel(sprintf('Ground Truth: r=%1.2f',r));
    
    subplot(length(Ins),3,(n-1)*3+3);
    Cn = reshape(Tnets(:,:,4,n),params.S*An,1); r = corr(Cg,Cn);
    scatplot(Cg,Cn); axis([-1 1 -1 1]); hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else plot(-1:0.1:1,1:-0.1:-1,'r'); end
    if n == 1; title(sprintf('Tnets scatterplot: %s (%s)',AllNames{4},ICA_Tnet),'color','r','interpreter','none'); end
    xlabel(sprintf('Ground Truth: r=%1.2f',r));
end
print(gcf,'-dpng','-r300','Results/Comparison_temporal_Scatter.png')

% Plot spatial netmat scatter plots 
Scorrs = mean(Scorrs,1); Scorrs = squeeze(mean(Scorrs,2));
figure; set(gcf,'Position',[0 0 1900 955],'PaperPositionMode','auto')
for n = 1:length(Ins)
    Cg = reshape(Snets(:,:,1,n),params.S*An,1);
    
    subplot(length(Ins),3,(n-1)*3+1);
    Cn = reshape(Snets(:,:,2,n),params.S*An,1); r = corr(Cg,Cn);
    scatplot(Cg,Cn); axis([-1 1 -1 1]); hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else plot(-1:0.1:1,1:-0.1:-1,'r'); end
    if n == 1; title(sprintf('Snets scatterplot: %s',AllNames{2}),'color','b','interpreter','none'); end
    xlabel(sprintf('Snet r=%1.2f\n Subject map r=%1.2f',r,Scorrs(2,n)));
    ylabel(sprintf('%s',Ins{n}),'interpreter','none')
    
    subplot(length(Ins),3,(n-1)*3+2);
    Cn = reshape(Snets(:,:,3,n),params.S*An,1); r = corr(Cg,Cn);
    scatplot(Cg,Cn); axis([-1 1 -1 1]); hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else plot(-1:0.1:1,1:-0.1:-1,'r'); end
    if n == 1; title(sprintf('Snets scatterplot: %s',AllNames{3}),'color',[0.6 0 1],'interpreter','none'); end
    xlabel(sprintf('Snet r=%1.2f\n Subject map r=%1.2f',r,Scorrs(3,n)));
    
    subplot(length(Ins),3,(n-1)*3+3);
    Cn = reshape(Snets(:,:,4,n),params.S*An,1); r = corr(Cg,Cn);
    scatplot(Cg,Cn); axis([-1 1 -1 1]); hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else plot(-1:0.1:1,1:-0.1:-1,'r'); end
    if n == 1; title(sprintf('Snets scatterplot: %s (%s)',AllNames{4},ICA_Snet),'color','r','interpreter','none'); end
    xlabel(sprintf('Snet r=%1.2f\n Subject map r=%1.2f',r,Scorrs(4,n)));
end
print(gcf,'-dpng','-r300','Results/Comparison_spatial_Scatter.png')

warning on

