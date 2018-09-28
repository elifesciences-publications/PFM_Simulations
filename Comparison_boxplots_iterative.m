clear all; close all; clc


% Only works if simulations have the same number of voxels, subjects, and modes

% Inputs to compare
Ins = {'LowO_MedM_MedT'};
save_figs = 1;

% Set paths
warning off
addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/
addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
addpath /vols/Scratch/janineb/HCP/DMN1200/PFM_Simulations/Overlap_functions/

% Intitialize variables
load(sprintf('Results/PFMsims_atlas_%s_%02d.mat',Ins{1},1),'atlasParams','params');
Ai = ones(params.iN); Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);
Signs = zeros(params.iN,10,5,length(Ins));
Orders = zeros(params.iN,10,5,length(Ins));
GroupMaps = zeros(atlasParams.V,params.iN,10,5,length(Ins));
Snets = zeros(params.S,10,4,length(Ins));
Tnets = zeros(params.S,10,8,length(Ins));
Scorrs = zeros(params.S,params.iN,10,4,length(Ins));
Tcorrs = zeros(params.S,params.iN,10,8,length(Ins));
Tnames = {'PFM clean TS','PFM PROFUMO nets','PFM DR TS',...
    'tICA DR','tICA cut','mICA cut','sICA DR','sICA DR thres'};
Mnames = {'PFM subject','tICA DR','mICA SR','sICA DR'};

%%% Loop over all comparison datasets and extract relevant information:
for n = 1:length(Ins)
    for d = 1:10
        fprintf('Loading data for input %d, itteration %d (%s)\n',n,d,Ins{n})
        
        % Load data
        load(sprintf('Results/PFMsims_atlas_%s_%02d.mat',Ins{n},d),...
            'Pg','P','A','D','atlasParams','params',...
            'tica_P1_DRnew','tica_A1_DRnew','ticaP1','ticaA1','ticaP1_DR','ticaA1_DR',...
            'micaP1','micaA1','mica_P1_DRnew');
        
        % Rerun DR
        [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres] = Melodic_DR(Ins{n},D,atlasParams,params);
        
        % Rerun PFM gathering
        [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(Ins{n},params);
        
        % Get all group maps and fix order & sign to match ground truth:
        pfmPg_new = pfmPg1_new;
        [C12,munkres_assign] = spatialcorr(pfmPg1_new,Pg);
        [i_pfm_new,~] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
        pfmPg_new = pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
        GroupMaps(:,:,d,2,n) = pfmPg_new; Signs(:,d,2,n) = sign_pfm_new; Orders(:,d,2,n) = i_pfm_new;
        
        ticaPg = ticaP1{1};
        [C12,munkres_assign] = spatialcorr(ticaPg,Pg);
        [i_tica,~] = find(munkres_assign==1); sign_tica = sign(C12(munkres_assign==1));
        ticaPg = ticaPg(:,i_tica).*repmat(sign_tica',atlasParams.V,1);
        GroupMaps(:,:,3,d,n) = ticaPg; Signs(:,d,3,n) = sign_tica; Orders(:,d,3,n) = i_tica;
        
        micaPg = micaP1{1};
        [C12,munkres_assign] = spatialcorr(micaPg,Pg);
        [i_ica_new,~] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
        micaPg = micaPg(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1);
        GroupMaps(:,:,4,d,n) = micaPg; Signs(:,d,4,n) = sign_ica_new; Orders(:,d,4,n) = i_ica_new;
        
        [C12,munkres_assign] = spatialcorr(sicaPg_new,Pg);
        [i_ica_new,~] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
        sicaPg_new = sicaPg_new(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1);
        GroupMaps(:,:,5,d,n) = sicaPg_new; Signs(:,d,5,n) = sign_ica_new; Orders(:,d,5,n) = i_ica_new;
        
        GroupMaps(:,:,1,d,n) = Pg;
        clear pfmPg_new C12 munkres_assign i_pfm_new sign_pfm_new ticaPg i_tica sign_tica i_ica_new sign_ica_new
        
        % Get subject temporal and spatial netmats
        Nedge = params.iN*params.iN;
        for s = 1:params.S
            
            % Ground truth nets
            [r,~] = spatialcorr(P{s},P{s}); SG = r(Ai);
            r = corr([A{s}{1}'; A{s}{2}']); TG = r(Ai);
            TSG = [A{s}{1}'; A{s}{2}'];
            
            % PFM nets
            p = pfmP1_new{s}; [r,~] = spatialcorr(p(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',atlasParams.V,1), p(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',atlasParams.V,1));
            Snets(s,d,1,n) = corr(r(Ai),SG);
            [r,~] = spatialcorr(p(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',atlasParams.V,1),P{s});
            Scorrs(s,:,d,1,n) = r(eye(params.iN)==1)';
            % PFM clean ts
            a = [pfmA1_new{s}{1}'; pfmA1_new{s}{2}'];
            r = corr(a(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,1,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',params.T*2,1));
            Tnets(s,d,1,n) = corr(r(Ai),TG);
            clear a
            % PFM PROFUMO
            [~,a] = cov2corr(inv(pfmNet{s}));
            a = a(:,Orders(:,d,2,n)); a = a(Orders(:,d,2,n),:);
            a = a.*(repmat(Signs(:,d,2,n)',params.iN,1) .* repmat(Signs(:,d,2,n),1,params.iN));
            Tnets(s,d,2,n) = corr(a(Ai),TG);
            clear a
            % PFM DR TS
            a = [];
            for r = 1:params.R(1)
                a = [a; nets_demean((pinv(nets_demean(p)) * nets_demean(double(D{s}{r})))')];
            end
            r = corr(a(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,3,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,2,n)).*repmat(Signs(:,d,2,n)',params.T*2,1));
            Tnets(s,d,3,n) = corr(r(Ai),TG);
            clear p a
            
            % tICA nets
            p = tica_P1_DRnew{s}; [r,~] = spatialcorr(p(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',atlasParams.V,1), p(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',atlasParams.V,1));
            Snets(s,d,2,n) = corr(r(Ai),SG);
            [r,~] = spatialcorr(p(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',atlasParams.V,1),P{s});
            Scorrs(s,:,d,2,n) = r(eye(params.iN)==1); clear p
            a = [tica_A1_DRnew{s}{1}; tica_A1_DRnew{s}{2}];
            r = corr(a(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,4,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',params.T*2,1));
            Tnets(s,d,4,n) = corr(r(Ai),TG); clear a
            
            a = [ticaA1{s}{1}'; ticaA1{s}{2}'];
            r = corr(a(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,5,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,3,n)).*repmat(Signs(:,d,3,n)',params.T*2,1));
            Tnets(s,d,5,n) = corr(r(Ai),TG); clear a
            
            % mICA nets
            p = mica_P1_DRnew{s}; [r,~] = spatialcorr(p(:,Orders(:,d,4,n)).*repmat(Signs(:,d,4,n)',atlasParams.V,1), p(:,Orders(:,d,4,n)).*repmat(Signs(:,d,4,n)',atlasParams.V,1));
            Snets(s,d,3,n) = corr(r(Ai),SG);
            [r,~] = spatialcorr(p(:,Orders(:,d,4,n)).*repmat(Signs(:,d,4,n)',atlasParams.V,1),P{s});
            Scorrs(s,:,d,3,n) = r(eye(params.iN)==1); clear p
            a = [micaA1{s}{1}'; micaA1{s}{2}'];
            r = corr(a(:,Orders(:,d,4,n)).*repmat(Signs(:,d,4,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,6,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,4,n)).*repmat(Signs(:,d,4,n)',params.T*2,1));
            Tnets(s,d,6,n) = corr(r(Ai),TG); clear a
            
            % sICA nets
            p = sica_P1_DR_thres{s};
            % sICA original spatial measures
            [r,~] = spatialcorr(p(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',atlasParams.V,1), p(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',atlasParams.V,1));
            Snets(s,d,4,n) = corr(r(Ai),SG);
            [r,~] = spatialcorr(p(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',atlasParams.V,1),P{s});
            Scorrs(s,:,d,4,n) = r(eye(params.iN)==1);
            % sICA original temporal measures
            a = [sica_A1_DR_new{s}{1}; sica_A1_DR_new{s}{2}];
            r = corr(a(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,7,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',params.T*2,1));
            Tnets(s,d,7,n) = corr(r(Ai),TG); clear a
            % sICA thresholded temporal measures
            a = [sica_A1_DR_thres{s}{1}; sica_A1_DR_thres{s}{2}];
            r = corr(a(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',params.T*2,1),TSG);
            Tcorrs(s,:,d,8,n) = r(eye(params.iN)==1);
            r = corr(a(:,Orders(:,d,5,n)).*repmat(Signs(:,d,5,n)',params.T*2,1));
            Tnets(s,d,8,n) = corr(r(Ai),TG); clear a
        end
        clear('pfmPg1_new','Pg','ticaP1','sicaPg_new','P','A',...
            'pfmP1_new','pfmA1_new','ticaP1_DR','ticaA1_DR',...
            'sica_P1_DR_new','sica_A1_DR_new','D','P');
    end
end

% General figure settings
Tcmap = zeros(7,3);
Tcmap(1:3,:) = repmat([0 0 1],3,1); Tcmap(4:5,:) = repmat([1 0 0],2,1); Tcmap(6,:) = [0 1 0]; Tcmap(7:8,:) = repmat([0.6 0 1],2,1);
Mcmap = [0 0 1; 1 0 0; 0 1 0; 0.6 0 1];
if length(Ins) == 1; TN = Tnames; MN = Mnames;
elseif length(Ins) == 2; TN = [Tnames Tnames]; MN = [Mnames Mnames]; Tcmap = [Tcmap; Tcmap];  Mcmap = [Mcmap; Mcmap];
elseif length(Ins) == 3; TN = [Tnames Tnames Tnames]; MN = [Mnames Mnames Mnames]; Tcmap = [Tcmap; Tcmap; Tcmap]; Mcmap = [Mcmap; Mcmap; Mcmap];   
elseif length(Ins) == 4; TN = [Tnames Tnames Tnames Tnames]; MN = [Mnames Mnames Mnames Mnames]; Tcmap = [Tcmap; Tcmap; Tcmap; Tcmap]; Mcmap = [Mcmap; Mcmap; Mcmap; Mcmap];
end
    
% Plot Full netmat figure
D = reshape(Tnets,params.S*10,8*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',TN,'color',Tcmap); XYrotalabel(60,0)
text(1:8:8*length(Ins),repmat(max(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
%vline(8.5:8:8*length(Ins),'k'); hline([max(median(D)) 0],':k')
title('Correlation of full netmats (temporal)');
if save_figs==1; print(gcf,'-dpng','-r300','Results/itt3_Boxplots_Tnets.png'); end

D = reshape(Snets,params.S*10,4*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',MN,'color',Mcmap); XYrotalabel(60,0)
text(0.5:4:4*length(Ins),repmat(max(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
%vline(4.5:4:4*length(Ins),'k'); hline([max(median(D)) 0],':k')
title('Correlation of spatial correlation matrices');
if save_figs==1; print(gcf,'-dpng','-r300','Results/itt4_Boxplots_Snets.png'); end

D = reshape(Scorrs,params.S*params.iN*10,4*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',MN,'color',Mcmap); XYrotalabel(60,0)
text(0.5:4:4*length(Ins),repmat(min(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
%vline(4.5:4:4*length(Ins),'k'); hline([max(median(D)) 0],':k')
title('Correlation of spatial maps');
if save_figs==1; print(gcf,'-dpng','-r300','Results/itt1_Boxplots_Scorrs.png'); end

D = reshape(Tcorrs,params.S*params.iN*10,8*length(Ins));
figure; set(gcf,'Position',[0 0 1800 500],'PaperPositionMode','auto')
subplot('Position',[0.07 0.25 0.86 0.65])
boxplot(D,'Labels',TN,'color',Tcmap); XYrotalabel(60,0)
text(1:8:8*length(Ins),repmat(min(D(:)),1,length(Ins)),Ins,'interpreter','none','fontsize',14)
%vline(8.5:8:8*length(Ins),'k'); hline([max(median(D)) 0],':k')
title('Correlation of timecourses');
if save_figs==1; print(gcf,'-dpng','-r300','Results/itt2_Boxplots_Tcorrs.png'); end

