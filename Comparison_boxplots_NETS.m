close all

% PFM maps missed (should ignore):
PFMmis = find(std(GroupMaps(:,:,2))<0.001);
Ai = ones(params.N); Ai(PFMmis,:) = 0; Ai(:,PFMmis) = 0;
Ai = triu(Ai,1); Ai = find(Ai); An = length(Ai);

% Intialise parameters
TG = zeros(params.S,length(Ai));
TI = zeros(params.S,length(Ai));
TIold = zeros(params.S,length(Ai));
TP = zeros(params.S,length(Ai));
SG = zeros(params.S,length(Ai));
SP = zeros(params.S,length(Ai));
SI = zeros(params.S,length(Ai));
SIold = zeros(params.S,length(Ai));

for s = 1:params.S
        
        % Ground truth temporal nets
        r = corr([A{s}{1}'; A{s}{2}']); TG(s,:) = r(Ai); clear r
        
        % DR original temporal nets
        a = [sica_A1_DR_new{s}{1}; sica_A1_DR_new{s}{2}];
        a = a(:,Orders(:,5,n)).*repmat(Signs(:,5,n)',params.T*2,1);
        r = corr(a); TIold(s,:) = r(Ai); clear r a
        
        % DR threshold temporal nets
        a = [sica_A1_DR_thres{s}{1}; sica_A1_DR_thres{s}{2}]; 
        a = a(:,Orders(:,5,n)).*repmat(Signs(:,5,n)',params.T*2,1);
        r = corr(a); TI(s,:) = r(Ai); clear r a
        
        % Ground truth spatial netmats
        [r,~] = spatialcorr(P{s},P{s}); SG(s,:) = r(Ai); clear r
        
        % PFM spatial nets
        p = pfmP1_new{s}; [r,~] = spatialcorr(p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1), p(:,Orders(:,2,n)).*repmat(Signs(:,2,n)',atlasParams.V,1)); 
        SP(s,:) = r(Ai); clear r p
        
        % PFM temporal nets
        [~,a] = cov2corr(inv(pfmNet{s}));
        a = a(:,Orders(:,2,n)); a = a(Orders(:,2,n),:);
        a = a.*(repmat(Signs(:,2,n)',params.N,1) .* repmat(Signs(:,2,n),1,params.N));
        TP(s,:) = a(Ai); clear a
        
        % ICA thresholded spatial nets
        p = sica_P1_DR_thres{s};
        [r,~] = spatialcorr(p(:,Orders(:,5,n)).*repmat(Signs(:,5,n)',atlasParams.V,1), p(:,Orders(:,5,n)).*repmat(Signs(:,5,n)',atlasParams.V,1));
        SI(s,:) = r(Ai); clear p r
        
        % ICA thresholded spatial nets
        p = sica_P1_DR_new{s};
        [r,~] = spatialcorr(p(:,Orders(:,5,n)).*repmat(Signs(:,5,n)',atlasParams.V,1), p(:,Orders(:,5,n)).*repmat(Signs(:,5,n)',atlasParams.V,1));
        SIold(s,:) = r(Ai); clear p r
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = 2;
[H,p,CI,Tstat] = ttest(SG);
EdgeSelectSpat = find(p<0.05/size(TI,2)); spat = [];

% Plot spatial results
for n = EdgeSelectSpat
    spat = [spat SG(:,n) SI(:,n) SIold(:,n) SP(:,n)];
end
figure; set(gcf,'Position',[0 0 1945 600],'PaperPositionMode','auto')
subplot('Position',[0.05 0.05 0.94 0.9])
boxplot(spat); hold on; vline(4.5:4:length(EdgeSelectSpat)*4,'k'); hline(0,'k');
scatter(1:4*length(EdgeSelectSpat),spat(s,:),15,'k','filled');
set(gca,'xtick',1:4*length(EdgeSelectSpat),'xticklabel',{'GT','ICAt','ICAo','PFM'})
title('Comparison of full spatial edges')
print(gcf,'-dpng','-r300','Results/Comparison_NETS_01.png')

spat = []; 
for n = EdgeSelectSpat
    [x,y] = ind2sub([params.iN params.iN],Ai(n));
    map_temp = pfmP1_new{s}(:,Orders(:,2,1)).*repmat(Signs(:,2,1)',atlasParams.V,1);
    spat = [spat P{s}(:,[x y]) 30*map_temp(:,[x y])];
end   
figure; set(gcf,'Position',[0 0 1945 600],'PaperPositionMode','auto')
subplot('Position',[0.02 0.05 0.98 0.9])
imagesc(spat,[-130 130]); hold on; vline(4.5:4:length(EdgeSelectSpat)*4,'k');
set(gca,'xtick',1.5:2:4*length(EdgeSelectSpat),'xticklabel',{'GT','PFM'})
title('Example subject ground truth and PFM maps')
print(gcf,'-dpng','-r300','Results/Comparison_NETS_03.png')

spat = [];
for n = EdgeSelectSpat
    [x,y] = ind2sub([params.iN params.iN],Ai(n));
    GM_thres = zeros(params.V,2,params.S); GM_dr = zeros(params.V,2,params.S);
    map_temp = sica_P1_DR_thres{s}(:,Orders(:,5,1)).*repmat(Signs(:,5,1)',atlasParams.V,1);
    spat = [spat map_temp(:,[x y])];
    map_temp = sica_P1_DR_new{s}(:,Orders(:,5,1)).*repmat(Signs(:,5,1)',atlasParams.V,1);
    spat = [spat map_temp(:,[x y])];
end
figure; set(gcf,'Position',[0 0 1945 600],'PaperPositionMode','auto')
subplot('Position',[0.02 0.05 0.98 0.9])
imagesc(spat,[-30 30]); hold on; vline(4.5:4:length(EdgeSelectSpat)*4,'k');
set(gca,'xtick',1.5:2:4*length(EdgeSelectSpat),'xticklabel',{'ICAthres','ICAdr'})
title('Example subject ICA DR maps before and after thresholding')
print(gcf,'-dpng','-r300','Results/Comparison_NETS_04.png')

% Plot temporal results
temp = []; 
for n = EdgeSelectSpat
    temp = [temp TG(:,n) TI(:,n) TIold(:,n) TP(:,n)];
end
figure; set(gcf,'Position',[0 0 1945 600],'PaperPositionMode','auto')
subplot('Position',[0.05 0.05 0.94 0.9])
boxplot(temp); hold on; vline(4.5:4:length(EdgeSelectSpat)*4,'k'); hline(0,'k');
scatter(1:4*length(EdgeSelectSpat),temp(s,:),15,'k','filled');
set(gca,'xtick',1:4*length(EdgeSelectSpat),'xticklabel',{'GT','ICAt','ICAo','PFM'})
title('Comparison of full temporal edges');
print(gcf,'-dpng','-r300','Results/Comparison_NETS_02.png')



