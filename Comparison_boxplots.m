clear all; close all; clc


% Only works if simulations have the same number of voxels, subjects, and modes

% Inputs to compare
Ins = {'original','HighO_LowM_LowT','LowO_LowM_HighT','LowO_HighM_LowT'};
    
% Set paths
warning off
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
addpath /vols/Scratch/janineb/HCP/DMN1200/PFM_Simulations/Overlap_functions/

% Intitialize variables
load(sprintf('Results/PFMsims_atlas_%s.mat',Ins{1}),'atlasParams','params');
Ai = ones(params.iN); Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);
Signs = zeros(params.iN,4,length(Ins));
Orders = zeros(params.iN,4,length(Ins));
GroupMaps = zeros(atlasParams.V,params.iN,4,length(Ins));
Snets = zeros(params.S,7,length(Ins));
Tnets = zeros(params.S,7,length(Ins));
Scorrs = zeros(params.S,params.iN,7,length(Ins));
Tcorrs = zeros(params.S,params.iN,7,length(Ins));
AllNames = {'PFM clean TS','PFM PROFUMO','PFM DR TS','tICA lowdim PCA','tICA highdim ICA','sICA group DR','sICA group DR thres'};
        
%%% Loop over all comparison datasets and extract relevant information:
for n = 1:length(Ins)
    fprintf('Loading data for input %d (%s)\n',n,Ins{n})
    
    % Load data
    load(sprintf('Results/PFMsims_atlas_%s.mat',Ins{n}),...
        'pfmPg1_new','Pg','ticaP1','P','A','sticaP1_DR','sticaA1_DR',...
        'pfmP1_new','pfmA1_new','ticaP1_DR','ticaA1_DR',...
        'D','atlasParams','params','P'); 
    
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
        [r,~] = spatialcorr(P{s},P{s}); SG = r(Ai);
        r = corr([A{s}{1}'; A{s}{2}']); TG = r(Ai);
        TSG = [A{s}{1}'; A{s}{2}'];
        
        % PFM nets
        p = pfmP1_new{s}; [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1), p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1)); 
        Snets(s,1,n) = corr(r(Ai),SG);
        [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1),P{s}); 
        Scorrs(s,:,1,n) = r(eye(params.iN)==1)'; 
        % PFM clean ts
            a = [pfmA1_new{s}{1}'; pfmA1_new{s}{2}']; 
            r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1),TSG); 
            Tcorrs(s,:,1,n) = r(eye(params.iN)==1); 
            r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1)); 
            Tnets(s,1,n) = corr(r(Ai),TG); 
            clear a
        % PFM PROFUMO
            [~,a] = cov2corr(inv(pfmNet{s})); 
            a = a(:,Orders(:,2,n)); a = a(Orders(:,2,n),:);
            a = a.*(repmat(Signs(:,2,n)',params.iN,1) .* repmat(Signs(:,2,n),1,params.iN));
            Tnets(s,2,n) = corr(a(Ai),TG); 
            clear a
        % PFM DR TS
            a = [];
            for r = 1:params.R(1)
                a = [a; nets_demean((pinv(nets_demean(p)) * nets_demean(double(D{s}{r})))')];
            end
            r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1),TSG); 
            Tcorrs(s,:,3,n) = r(eye(params.iN)==1); 
            r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1)); 
            Tnets(s,3,n) = corr(r(Ai),TG); 
            clear p a
        
        % tICA nets
        p = ticaP1_DR{s}; [r,~] = spatialcorr(p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1), p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1)); 
        Snets(s,4,n) = corr(r(Ai),SG);
        [r,~] = spatialcorr(p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1),P{s}); 
        Scorrs(s,:,4,n) = r(eye(params.iN)==1); clear p
        a = [ticaA1_DR{s}{1}'; ticaA1_DR{s}{2}']; 
        r = corr(a(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',params.T*2,1),TSG); 
        Tcorrs(s,:,4,n) = r(eye(params.iN)==1); 
        r = corr(a(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',params.T*2,1)); 
        Tnets(s,4,n) = corr(r(Ai),TG); clear a
        
        p = sticaP1_DR{s}; [r,~] = spatialcorr(p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1), p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1)); 
        Snets(s,5,n) = corr(r(Ai),SG);
        [r,~] = spatialcorr(p(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',atlasParams.V,1),P{s}); 
        Scorrs(s,:,5,n) = r(eye(params.iN)==1); clear p
        a = [sticaA1_DR{s}{1}'; sticaA1_DR{s}{2}']; 
        r = corr(a(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',params.T*2,1),TSG); 
        Tcorrs(s,:,5,n) = r(eye(params.iN)==1); 
        r = corr(a(:,Orders(:,3,n)).*repmat(Signs(:,3,n)',params.T*2,1)); 
        Tnets(s,5,n) = corr(r(Ai),TG); clear a
        
        % sICA nets
        p = sica_P1_DR_thres{s}; 
        % sICA original spatial measures
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1), p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1)); 
            Snets(s,6,n) = corr(r(Ai),SG); 
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1),P{s}); 
            Scorrs(s,:,6,n) = r(eye(params.iN)==1);
        % sICA thresholded spatial measures
            p(abs(p)<12) = 0; 
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1), p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1)); 
            Snets(s,7,n) = corr(r(Ai),SG); 
            [r,~] = spatialcorr(p(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',atlasParams.V,1),P{s}); 
            Scorrs(s,:,7,n) = r(eye(params.iN)==1); clear p
        % sICA original temporal measures         
            a = [sica_A1_DR_new{s}{1}; sica_A1_DR_new{s}{2}]; 
            r = corr(a(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',params.T*2,1),TSG); 
            Tcorrs(s,:,6,n) = r(eye(params.iN)==1); 
            r = corr(a(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',params.T*2,1)); 
            Tnets(s,6,n) = corr(r(Ai),TG); clear a
        % sICA thresholded temporal measures         
            a = [sica_A1_DR_thres{s}{1}; sica_A1_DR_thres{s}{2}]; 
            r = corr(a(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',params.T*2,1),TSG); 
            Tcorrs(s,:,7,n) = r(eye(params.iN)==1); 
            r = corr(a(:,Orders(:,4,n)).*repmat(Signs(:,4,n)',params.T*2,1)); 
            Tnets(s,7,n) = corr(r(Ai),TG); clear a
    end
    clear('pfmPg1_new','Pg','ticaP1','sicaPg_new','P','A',...
        'pfmP1_new','pfmA1_new','ticaP1_DR','ticaA1_DR',...
        'sica_P1_DR_new','sica_A1_DR_new','D','P'); 
end

% General figure settings
cmap = zeros(7,3);
cmap(1:3,:) = repmat([0 0 1],3,1); cmap(4:5,:) = repmat([1 0 0],2,1); cmap(6:7,:) = repmat([0.6 0 1],2,1);
if length(Ins) == 1; N = AllNames;
elseif length(Ins) == 2; N = [AllNames AllNames]; cmap = [cmap; cmap];
elseif length(Ins) == 3; N = [AllNames AllNames AllNames]; cmap = [cmap; cmap; cmap];   
elseif length(Ins) == 4; N = [AllNames AllNames AllNames AllNames]; cmap = [cmap; cmap; cmap; cmap];
end
    
% Plot Full netmat figure
D = reshape(Tnets,params.S,7*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',N,'color',cmap); XYrotalabel(60,0)
text(3.05:7:7*length(Ins),repmat(max(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
vline(7.5:7:7*length(Ins),'k'); hline(max(median(D)),':k')
title('Correlation of full netmats (temporal)');
print(gcf,'-dpng','-r300','Results/Boxplots_Tnets.png')

D = reshape(Snets,params.S,7*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',N,'color',cmap); XYrotalabel(60,0)
text(3.05:7:7*length(Ins),repmat(max(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
vline(7.5:7:7*length(Ins),'k'); hline(max(median(D)),':k')
title('Correlation of spatial correlation matrices');
print(gcf,'-dpng','-r300','Results/Boxplots_Snets.png')

D = reshape(Scorrs,params.S*params.iN,7*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',N,'color',cmap); XYrotalabel(60,0)
text(3.05:7:7*length(Ins),repmat(min(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
vline(7.5:7:7*length(Ins),'k'); hline(max(median(D)),':k')
title('Correlation of spatial maps');
print(gcf,'-dpng','-r300','Results/Boxplots_Scorrs.png')

D = reshape(Tcorrs,params.S*params.iN,7*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',N,'color',cmap); XYrotalabel(60,0)
text(3.05:7:7*length(Ins),repmat(min(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
vline(7.5:7:7*length(Ins),'k'); hline(max(median(D)),':k')
title('Correlation of timecourses');
print(gcf,'-dpng','-r300','Results/Boxplots_Tcorrs.png')

Ns = [3 4 7];
figure; set(gcf,'Position',[0 0 1500 955],'PaperPositionMode','auto'); I = 1;
c = [0 0 1; 1 0 0; 0.6 0 1];
for m = 1:4
    for n = 1:3
        subplot(4,3,I) 
        boxplot(squeeze(Tcorrs(:,:,Ns(n),m))); I = I+1;
        axis([0.5 params.iN+0.5 -0.6 1.1])
        if m == 1; title(sprintf('timecourses %s',AllNames{Ns(n)}),'interpreter','none','color',c(n,:)); end
        if n == 1; ylabel(sprintf('%s',Ins{m}),'interpreter','none'); end
    end
end
print(gcf,'-dpng','-r300','Results/Boxplots_Tcorrs_separate.png')

Ns = [1 4 7];
figure; set(gcf,'Position',[0 0 1500 955],'PaperPositionMode','auto'); I = 1;
c = [0 0 1; 1 0 0; 0.6 0 1];
for m = 1:4
    for n = 1:3
        subplot(4,3,I) 
        boxplot(squeeze(Scorrs(:,:,Ns(n),m))); I = I+1;
        axis([0.5 params.iN+0.5 -0.6 1.1])
        if m == 1; title(sprintf('spatial maps %s',AllNames{Ns(n)}),'interpreter','none','color',c(n,:)); end
        if n == 1; ylabel(sprintf('%s',Ins{m}),'interpreter','none'); end
    end
end
print(gcf,'-dpng','-r300','Results/Boxplots_Scorrs_separate.png')
