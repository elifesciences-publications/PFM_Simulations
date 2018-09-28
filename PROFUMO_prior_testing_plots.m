clear all; close all; clc

% Rerun PFM gathering
Ins = {'MedO_MedM_MedT_01','MedO_MedM_MedT_01_20','MedO_MedM_MedT_01_strong20','MedO_MedM_MedT_01_med20','MedO_MedM_MedT_01_weak20'};

% Load data
load(sprintf('Results/PFMsims_atlas_%s.mat',Ins{1}),...
    'Pg','P','A','D','atlasParams','params');

% Set outputs
Ai = ones(params.N); Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);
Signs = zeros(params.N,5,length(Ins));
Orders = zeros(params.N,5,length(Ins));
GroupMaps = zeros(atlasParams.V,params.N,length(Ins)+1);
TEdges = zeros(params.S,length(Ai),length(Ins)+1);
SEdges = zeros(params.S,length(Ai),length(Ins)+1);
Snets = zeros(params.S,4,length(Ins));
Tnets = zeros(params.S,8,length(Ins));
Scorrs = zeros(params.S,params.N,4,length(Ins));
Tcorrs = zeros(params.S,params.N,8,length(Ins));
Tnames = {'PFM clean TS','PFM PROFUMO nets','PFM DR TS',...
    'tICA DR','tICA cut','mICA cut','sICA DR','sICA DR thres'};
Mnames = {'PFM subject','tICA DR','mICA SR','sICA DR'};

% GT group maps
GroupMaps(:,:,length(Ins)+1) = Pg;

for n = 1:length(Ins)
    [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(Ins{n},params);
    
    pfmPg_new = pfmPg1_new;
    [C12,munkres_assign] = spatialcorr(pfmPg1_new,Pg);
    [i_pfm_new,~] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
    pfmPg_new = pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1);
    GroupMaps(:,:,n) = pfmPg_new; Signs(:,2,n) = sign_pfm_new; Orders(:,2,n) = i_pfm_new;
    
    for s = 1:params.S
        
        % Ground truth nets
        [r,~] = spatialcorr(P{s},P{s}); 
        SEdges(s,:,5) = r(Ai); SG = r(Ai);
        r = corr([A{s}{1}'; A{s}{2}']); 
        TEdges(s,:,5) = r(Ai); TG = r(Ai);
        TSG = [A{s}{1}'; A{s}{2}'];
        
        % PFM nets
        p = pfmP1_new{s}; [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1), p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1));
        SEdges(s,:,n) = r(Ai);
        Snets(s,1,n) = corr(r(Ai),SG);
        [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1),P{s});
        Scorrs(s,:,1,n) = r(eye(params.N)==1)';
        % PFM clean ts
        a = [pfmA1_new{s}{1}'; pfmA1_new{s}{2}'];
        r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1),TSG);
        Tcorrs(s,:,1,n) = r(eye(params.N)==1);
        r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1));
        Tnets(s,1,n) = corr(r(Ai),TG);
        clear a
        % PFM PROFUMO
        [~,a] = cov2corr(inv(pfmNet{s}));
        a = a(:,Orders(:,2,n)); a = a(Orders(:,2,n),:);
        a = a.*(repmat(Signs(:,2,n)',params.N,1) .* repmat(Signs(:,2,n),1,params.N));
        TEdges(s,:,n) = a(Ai);
        Tnets(s,2,n) = corr(a(Ai),TG);
        clear a
        % PFM DR TS
        a = [];
        for r = 1:params.R(1)
            a = [a; nets_demean((pinv(nets_demean(p)) * nets_demean(double(D{s}{r})))')];
        end
        r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1),TSG);
        Tcorrs(s,:,3,n) = r(eye(params.N)==1);
        r = corr(a(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',params.T*2,1));
        Tnets(s,3,n) = corr(r(Ai),TG);
        clear p a
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

figure; set(gcf,'Position',[0 570 1600 500],'PaperPositionMode','auto')
for n = 1:length(Ins)+1
    subplot(1,length(Ins)+1,n); imagesc(GroupMaps(:,:,n)); colormap parula; set(gca,'xtick',1:params.N); grid on
    if n == length(Ins)+1; ylabel('Ground truth group maps'); else ylabel(sprintf('%s',Ins{n}),'interpreter','none'); end
end

figure; 
subplot(1,2,1); plot(squeeze(mean(TEdges,1)));
subplot(1,2,2); plot(squeeze(mean(SEdges,1)));
legend([Ins {'GT'}],'interpreter','none')

