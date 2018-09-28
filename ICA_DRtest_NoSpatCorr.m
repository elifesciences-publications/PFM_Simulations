clear all; close all; clc

addpath /Users/janineb/Downloads/

% Set parameters & initiate variables
V=10000; t=200; D = 2; S = 50;
Y = zeros(t*S,V);
GTmaps = zeros(V,D*S);
GTts = zeros(t,D*S);
GTcorrs = zeros(S,2);
DRcorrs = zeros(S,2);
MapTcorrs = zeros(S,4);

% Generate data
for s = 1:S
    A = laprnd(V,D,0,0.5);
    A(1000:1100,1) = A(1000:1100,1)+rand(101,1)*10+2;
    A(1050:1150,2) = A(1050:1150,2)+rand(101,1)*10+2;
    T = nets_demean(randn(t,D));
    T = T+repmat(randn(t,1),1,2); % Add temporal correlation
    GTcorrs(s,1) = corr(T(:,1),T(:,2)); GTcorrs(s,2) = corr(A(:,1),A(:,2));
    GTmaps(:,(s-1)*D+1:s*D) = A; GTts(:,(s-1)*D+1:s*D) = T;
    Y((s-1)*t+1:s*t,:) = nets_demean( T * A' + randn(t,V)*0.1 );
end

% Run group ICA
[icaS,icaA,icaW] = fastica(Y,'approach', 'symm', 'g', 'tanh', 'lastEig',2);

% Reorder and change sign of IC maps to match GT
R = corr([mean(reshape(GTmaps,V,D,S),3) icaS']); R = R(3:end,1:2);
[~,i] = max(abs(R)); if length(unique(i))==1; noi = setdiff(1:2,i); [~,inew] = max(abs(R(noi,:))); i(inew) = noi; end
icaS = icaS(i,:);
R = corr([mean(reshape(GTmaps,V,D,S),3) icaS']); R = R(3:end,1:2);
icaS(1,:) = icaS(1,:)*sign(R(1,1)); icaS(2,:) = icaS(2,:)*sign(R(2,2)); 

% Run dual regression
pGM = pinv(nets_demean(icaS'));
DRmaps = zeros(V,D,S);
for s = 1:S
    d = Y((s-1)*t+1:s*t,:)';
    NODEts = nets_demean((pGM*nets_demean(d))');
    MAP = ssglmT(nets_demean(d'),NODEts)';
    DRmaps(:,:,s) = MAP;
    DRcorrs(s,1) = corr(NODEts(:,1),NODEts(:,2)); DRcorrs(s,2) = corr(MAP(:,1),MAP(:,2)); 
    MapTcorrs(s,1) = corr(MAP(:,1),GTmaps(:,(s-1)*D+1));
    MapTcorrs(s,2) = corr(MAP(:,2),GTmaps(:,(s-1)*D+2));
    MapTcorrs(s,3) = corr(NODEts(:,1),GTts(:,(s-1)*D+1));
    MapTcorrs(s,4) = corr(NODEts(:,2),GTts(:,(s-1)*D+2));
end

% Save Results
save('ICA_DRtest_NoSpatCorr.mat');

% Plot Results
figure; set(gcf,'Position',[1150 450 800 600],'PaperPositionMode','auto')
M = [GTcorrs DRcorrs MapTcorrs];
imagesc(corr(M),[-1 1]); colorbar
hline([2.5 4.5],'k'); vline([2.5 4.5],'k');
N = {'GT Tnets','GT Snets','DR Tnets','DR Snets',...
    'Map1 corr','Map2 corr','TS1 corr','TS2 corr'};
set(gca,'xtick',1:8,'xticklabel',N,'ytick',1:8,'yticklabel',N,'FontSize',14)
for n = 1:size(M,2); text(7.7,n,sprintf('%1.4f',mean(M(:,n))),'FontSize',14); end
title('correlations across subjects (no GT spatial correlation)','FontSize',16);
print('SimulationDR_NoSpatCorr_maps','-dpng','-r150')

figure; set(gcf,'Position',[150 75 900 1000],'PaperPositionMode','auto')
M = [mean(reshape(GTmaps,V,D,S),3) icaS' mean(DRmaps,3)];
M = M.*repmat(sign(mean(M)),size(M,1),1)./repmat(max(abs(M)),size(M,1),1);
imagesc(M);
vline([2.5 4.5],'k'); colorbar
R = corr(mean(reshape(GTmaps,V,D,S),3)); text(1,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
R = corr(icaS'); text(3,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
R = corr(mean(DRmaps,3)); text(5,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
N = {'GT map1','GT map2','gICA map1','gICA map2','DR map1','DR map2'};
set(gca,'xtick',1:6,'xticklabel',N,'FontSize',14)
title('Group average maps (no GT spatial correlation)','FontSize',16)
print('SimulationDR_NoSpatCorr_comparisons','-dpng','-r150')

figure;
M = [GTcorrs DRcorrs MapTcorrs];
N = {'GT Tnets','GT Snets','DR Tnets','DR Snets',...
    'Map1 corr','Map2 corr','TS1 corr','TS2 corr'};
for n = 1:4
    subplot(2,2,n); scatter(M(:,6),M(:,n));
    xlabel(N{6}); ylabel(N{n});
end

figure;
[r,i] = sort(MapTcorrs(:,2));
m = reshape(GTmaps,V,D,S);
M = [m(:,2,i(1)) DRmaps(:,2,i(1)) m(:,2,i(end)) DRmaps(:,2,i(end))];
M = M.*repmat(sign(mean(M)),size(M,1),1)./repmat(max(abs(M)),size(M,1),1);
imagesc(M);


