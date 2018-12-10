clear all; close all; clc

Ins = {'New_sims_NoOverlap_NEWICA'}; %{'New_sims_HighM'};
warning off

% Set paths
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
addpath /vols/Scratch/janineb/HCP/DMN1200/PFM_Simulations/Overlap_functions/

% Load data
load(sprintf('Results/PFMsims_atlas_%s_01.mat',Ins{1}),'atlasParams','params');

% Intitialize variables
Ai = ones(params.N); Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);
rho = 0.1;
Signs = zeros(params.N,2,params.nRepeats);
Orders = zeros(params.N,2,params.nRepeats);
GroupMaps = zeros(atlasParams.V,params.N,3,params.nRepeats);
Snets = zeros(params.S,3,params.nRepeats);
Snets_subs = nan(length(Ai),3,params.nRepeats);
Tnets = zeros(params.S,4,params.nRepeats);
Tnets_subs = nan(length(Ai),4,params.nRepeats);
Scorrs = zeros(params.S,params.N,3,params.nRepeats);
Tcorrs = zeros(params.S,params.N,4,params.nRepeats);
Tnetmats = zeros(params.S,length(Ai),5,params.nRepeats);
Snetmats = zeros(params.S,length(Ai),4,params.nRepeats);

for d = 1:params.nRepeats
     
    filename = sprintf('%s_%02d',Ins{1},d);
    
    % Load data
    load(sprintf('Results/PFMsims_atlas_%s.mat',filename),...
        'Pg','P','A','D','atlasParams','params');
    fprintf('Running iteration %d out of %d\n',d,params.nRepeats);
    
    % Rerun DR
    [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres_median,sica_P1_DR_thres_null1,sica_P1_DR_thres_null2,sica_P1_DR_thres_intersect,sica_P1_DR_thres_null3,sica_A1_DR_thres_median,sica_A1_DR_thres_intersect,sica_A1_DR_thres_null1,sica_A1_DR_thres_null2,sica_A1_DR_thres_null3] = Melodic_DR(filename,D,atlasParams,params);
    clear D
    
    % Rerun PFM gathering
    [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(filename,params);
    
    % Get all group maps and fix order & sign to match ground truth:
    pfmPg_new = pfmPg1_new;
    [C12,munkres_assign] = spatialcorr(pfmPg1_new,Pg);
    [i_pfm_new,~] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
    pfmPg_new = pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
    GroupMaps(:,:,1,d) = pfmPg_new; Signs(:,1,d) = sign_pfm_new; Orders(:,1,d) = i_pfm_new;
    
    [C12,munkres_assign] = spatialcorr(sicaPg_new,Pg);
    [i_ica_new,~] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
    sicaPg_new = sicaPg_new(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1);
    GroupMaps(:,:,2,d) = sicaPg_new; Signs(:,2,d) = sign_ica_new; Orders(:,2,d) = i_ica_new;
    
    GroupMaps(:,:,3,d) = Pg;
    clear pfmPg_new C12 munkres_assign i_pfm_new sign_pfm_new ticaPg i_tica sign_tica i_ica_new sign_ica_new;
    
    % Get subject temporal and spatial netmats
    for s = 1:params.S
        
        % Ground truth nets
        [r,~] = corr(P{s},P{s}); Snetmats(s,:,4,d) = r(Ai); SG = r(Ai); 
        grot = cov([A{s}{1}'; A{s}{2}']);  
        grot = grot/sqrt(mean(diag(grot).^2)); grot = inv(grot+rho*eye(params.iN));
        grot = -grot; grot = (grot ./ repmat(sqrt(abs(diag(grot))),1,params.iN)) ./ repmat(sqrt(abs(diag(grot)))',params.iN,1); 
        Tnetmats(s,:,4,d) = grot(Ai); TG = grot(Ai);
        TSG = [A{s}{1}'; A{s}{2}'];
        clear r grot
        
        % PFM spatial nets
        p = pfmP1_new{s}; [r,~] = corr(p(:,Orders(:,1,d)).*repmat(Signs(:,1,d)',atlasParams.V,1), p(:,Orders(:,1,d)).*repmat(Signs(:,1,d)',atlasParams.V,1));
        Snets(s,1,d) = corr(r(Ai),SG); Snetmats(s,:,1,d) = r(Ai);
        r = corr(p(:,Orders(:,1,d)).*repmat(Signs(:,1,d)',atlasParams.V,1),P{s});
        Scorrs(s,:,1,d) = r(eye(params.N)==1)';
        clear p r 
        
        % PFM timeseries
        a = [pfmA1_new{s}{1}'; pfmA1_new{s}{2}'];
        r = corr(a(:,Orders(:,1,d)).*repmat(Signs(:,1,d)',params.T*2,1),TSG);
        Tcorrs(s,:,1,d) = r(eye(params.N)==1);
        
        % PFM PROFUMO
        [~,a] = cov2corr(inv(pfmNet{s}));
        a = a(:,Orders(:,1,d)); a = a(Orders(:,1,d),:);
        a = a.*(repmat(Signs(:,1,d)',params.N,1) .* repmat(Signs(:,1,d),1,params.N));
        Tnetmats(s,:,1,d) = a(Ai);
        Tnets(s,1,d) = corr(a(Ai),TG);
        clear a
        
        % sICA original spatial measures
        p = sica_P1_DR_new{s};
        [r,~] = corr(p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1), p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1));
        Snets(s,2,d) = corr(r(Ai),SG);
        Snetmats(s,:,2,d) = r(Ai);
        r = corr(p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1),P{s});
        Scorrs(s,:,2,d) = r(eye(params.N)==1);
        clear p r
        
        % sICA thresholded spatial measures (NULL)
        p = sica_P1_DR_thres_null2{s};
        r = corr(p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1), p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1));
        Snets(s,3,d) = corr(r(Ai),SG);
        Snetmats(s,:,3,d) = r(Ai);
        r = corr(p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1),P{s});
        Scorrs(s,:,3,d) = r(eye(params.N)==1);
        clear p r
        
        % sICA thresholded spatial measures (INTERSECT)
        p = sica_P1_DR_thres_null3{s};
        r = corr(p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1), p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1));
        Snets(s,4,d) = corr(r(Ai),SG);
        Snetmats(s,:,5,d) = r(Ai);
        r = corr(p(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',atlasParams.V,1),P{s});
        Scorrs(s,:,4,d) = r(eye(params.N)==1);
        clear p r
        
        % sICA original temporal measures
        a = [sica_A1_DR_new{s}{1}; sica_A1_DR_new{s}{2}];
        r = corr(a(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',params.T*2,1),TSG);
        Tcorrs(s,:,2,d) = r(eye(params.N)==1);
        grot = cov(a(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',params.T*2,1));  
        grot = grot/sqrt(mean(diag(grot).^2)); grot = inv(grot+rho*eye(params.iN));
        grot = -grot; grot = (grot ./ repmat(sqrt(abs(diag(grot))),1,params.iN)) ./ repmat(sqrt(abs(diag(grot)))',params.iN,1); 
        Tnetmats(s,:,2,d) = grot(Ai);
        Tnets(s,2,d) = corr(grot(Ai),TG); 
        clear a r grot
        
        % sICA thresholded temporal measures (NULL)
        a = [sica_A1_DR_thres_null2{s}{1}; sica_A1_DR_thres_null2{s}{2}];
        r = corr(a(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',params.T*2,1),TSG);
        Tcorrs(s,:,3,d) = r(eye(params.N)==1);
        grot = cov(a(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',params.T*2,1));  
        grot = grot/sqrt(mean(diag(grot).^2)); grot = inv(grot+rho*eye(params.iN));
        grot = -grot; grot = (grot ./ repmat(sqrt(abs(diag(grot))),1,params.iN)) ./ repmat(sqrt(abs(diag(grot)))',params.iN,1); 
        Tnetmats(s,:,3,d) = grot(Ai);
        Tnets(s,3,d) = corr(grot(Ai),TG); 
        clear a r grot
        
        % sICA thresholded temporal measures (INTERSECT)
        a = [sica_A1_DR_thres_null3{s}{1}; sica_A1_DR_thres_null3{s}{2}];
        r = corr(a(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',params.T*2,1),TSG);
        Tcorrs(s,:,4,d) = r(eye(params.N)==1);
        grot = cov(a(:,Orders(:,2,d)).*repmat(Signs(:,2,d)',params.T*2,1));  
        grot = grot/sqrt(mean(diag(grot).^2)); grot = inv(grot+rho*eye(params.iN));
        grot = -grot; grot = (grot ./ repmat(sqrt(abs(diag(grot))),1,params.iN)) ./ repmat(sqrt(abs(diag(grot)))',params.iN,1); 
        Tnetmats(s,:,5,d) = grot(Ai);
        Tnets(s,4,d) = corr(grot(Ai),TG); 
        clear a r grot
        
    end
    
    % Do edge correlations with GT across subjects (not across vectorised nets)
    
    for e = 1:length(Ai)
        [ix,iy] = ind2sub([params.iN params.iN],Ai(e));
        x = mean(mean(Scorrs(:,ix,:,d),1),3);
        y = mean(mean(Scorrs(:,iy,:,d),1),3);
        if min([squeeze(mean(Scorrs(:,ix,:,d),1)); squeeze(mean(Scorrs(:,iy,:,d),1))])>=0.5
            Tnets_subs(e,1,d) = corr(Tnetmats(:,e,1,d),Tnetmats(:,e,4,d));
            Tnets_subs(e,2,d) = corr(Tnetmats(:,e,2,d),Tnetmats(:,e,4,d));
            Tnets_subs(e,3,d) = corr(Tnetmats(:,e,3,d),Tnetmats(:,e,4,d));
            Tnets_subs(e,4,d) = corr(Tnetmats(:,e,5,d),Tnetmats(:,e,4,d));
            Snets_subs(e,1,d) = corr(Snetmats(:,e,1,d),Snetmats(:,e,4,d));
            Snets_subs(e,2,d) = corr(Snetmats(:,e,2,d),Snetmats(:,e,4,d));
            Snets_subs(e,3,d) = corr(Snetmats(:,e,3,d),Snetmats(:,e,4,d));
            Snets_subs(e,4,d) = corr(Snetmats(:,e,5,d),Snetmats(:,e,4,d));
        end
%         if min(squeeze(mean(Scorrs(:,ix,:,d),1)))<0.5
%             Tcorrs(:,ix,:,d) = nan;
%             Scorrs(:,ix,:,d) = nan;
%         end
%         if min(squeeze(mean(Scorrs(:,iy,:,d),1)))<0.5
%             Tcorrs(:,iy,:,d) = nan;
%             Scorrs(:,iy,:,d) = nan;
%         end
    end
    
    clear Pg P A D atlasParams 
    clear sicaPg_new sica_P1_DR_new sica_A1_DR_new sica_P1_DR_thres sica_A1_DR_thres
    clear pfmPg1_new pfmP1_new pfmA1_new pfmNet
end

% Save results and then plot
save(sprintf('Results/Paper_figure_results_%s.mat',Ins{1}),'Signs','Orders','GroupMaps',...
    'Snets','Tnets','Snets_subs','Tnets_subs','Scorrs','Tcorrs',...
    'Tnetmats','Snetmats','-v7.3');
warning on

figure; set(gcf,'Position',[0 0 850 945],'PaperPositionMode','auto')

% Plot correlations between timecourses
subplot(3,2,1);
histogram(Tcorrs(:,:,2,:),20,'DisplayStyle','stairs','edgecolor','g'); hold on;
histogram(Tcorrs(:,:,3,:),20,'DisplayStyle','stairs','edgecolor','b');
%histogram(Tcorrs(:,:,4,:),20,'DisplayStyle','stairs','edgecolor','k');
histogram(Tcorrs(:,:,1,:),20,'DisplayStyle','stairs','edgecolor','r');
legend({'Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'},'Location','NorthWest')
title('A. Accuracy of timeseries estimation');
xlabel(sprintf('Correlation between ground truth\n and estimated timeseries'))
axis([ -1 1 0 2000])

% Plot correlations between maps
subplot(3,2,2);
histogram(Scorrs(:,:,2,:),20,'DisplayStyle','stairs','edgecolor','g'); hold on;
histogram(Scorrs(:,:,3,:),20,'DisplayStyle','stairs','edgecolor','b');
%histogram(Scorrs(:,:,4,:),20,'DisplayStyle','stairs','edgecolor','k');
histogram(Scorrs(:,:,1,:),20,'DisplayStyle','stairs','edgecolor','r');
%legend({'Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'})
title('B. Accuracy of spatial map estimation');
xlabel(sprintf('Correlation between ground truth\n and estimated spatial map'))
axis([ -1 1 0 2000])

% Plot correlations between temporal edges and ground truth
subplot(3,2,3);
histogram(Tnets_subs(:,2,:),20,'DisplayStyle','stairs','edgecolor','g'); hold on;
histogram(Tnets_subs(:,3,:),20,'DisplayStyle','stairs','edgecolor','b');
%histogram(Tnets_subs(:,4,:),20,'DisplayStyle','stairs','edgecolor','k');
histogram(Tnets_subs(:,1,:),20,'DisplayStyle','stairs','edgecolor','r');
%legend({'Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'})
title('C. Accuracy of temporal edge estimation')
xlabel(sprintf('Correlation between ground truth\n and estimated temporal edges'))
axis([ -1 1 0 180])

% Plot correlations between spatial edges and ground truth
subplot(3,2,4);
histogram(Snets_subs(:,2,:),20,'DisplayStyle','stairs','edgecolor','g'); hold on;
histogram(Snets_subs(:,3,:),20,'DisplayStyle','stairs','edgecolor','b');
%histogram(Snets_subs(:,4,:),20,'DisplayStyle','stairs','edgecolor','k');
histogram(Snets_subs(:,1,:),20,'DisplayStyle','stairs','edgecolor','r');
%legend({'Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'})
title('D. Accuracy of spatial edge estimation')
xlabel(sprintf('Correlation between ground truth\n and estimated spatial edges'))
axis([ -1 1 0 180])

% Strong positive spatial edges only
[H,P,CI,STATS] = ttest(reshape(Snetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats),0,'tail','right');
I = find(P<0.05/size(P,2));
subplot(3,2,6);
A1 = reshape(Snetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A1 = A1(:,I);
A2 = reshape(Snetmats(:,:,2,:)-Snetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A2 = A2(:,I);
A3 = reshape(Snetmats(:,:,3,:)-Snetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A3 = A3(:,I);
A5 = reshape(Snetmats(:,:,5,:)-Snetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A5 = A5(:,I);
A4 = reshape(Snetmats(:,:,1,:)-Snetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A4 = A4(:,I);
boxplot([A1(:) A2(:) A3(:) A4(:)]); hold on; hline(0,':k'); vline(1.5,'k')
set(gca,'xtick',1:4,'xticklabel',{'Ground Truth','Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'})
title('F. Spatial edges (selected edges only)');
XYrotalabel(50); %ylabel('Edge estimate')
subplot(3,2,5);
A1 = reshape(Tnetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A1 = A1(:,I);
A2 = reshape(Tnetmats(:,:,2,:)-Tnetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A2 = A2(:,I);
A3 = reshape(Tnetmats(:,:,3,:)-Tnetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A3 = A3(:,I);
A5 = reshape(Tnetmats(:,:,5,:)-Tnetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A5 = A5(:,I);
A4 = reshape(Tnetmats(:,:,1,:)-Tnetmats(:,:,4,:),params.S,length(Ai)*params.nRepeats); A4 = A4(:,I);
boxplot([A1(:) A2(:) A3(:) A4(:)]); hold on; hline(0,':k'); vline(1.5,'k')
set(gca,'xtick',1:4,'xticklabel',{'Ground Truth','Traditional ICA-DR','Thresholded ICA-DR','PROFUMO'})
title('E. Temporal edges (selected edges only)');
 %ylabel('Edge estimate')

print(gcf,'-dpng','-r300',sprintf('Results/Paper_figure_%s.png',Ins{1}))



