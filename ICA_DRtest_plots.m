% Reorder and change sign of PFM maps to match GT
R = corr([mean(reshape(GTmaps,V,D,S),3) PFMgmap]); R = R(3:end,1:2);
[~,i] = max(abs(R)); if length(unique(i))==1; noi = setdiff(1:2,i); [~,inew] = max(abs(R(noi,:))); i(inew) = noi; end
PFMgmap = PFMgmap(:,i);
R = corr([mean(reshape(GTmaps,V,D,S),3) PFMgmap]); R = R(3:end,1:2);
PFMgmap(:,1) = PFMgmap(:,1)*sign(R(1,1)); PFMgmap(:,2) = PFMgmap(:,2)*sign(R(2,2));
for s = 1:S
    M = PFMmaps{s}; M = M(:,i);
    M(:,1) = M(:,1)*sign(R(1,1)); M(:,2) = M(:,2)*sign(R(2,2));
    PFMmaps{s} = M; clear M
    M = PFMts{s}; M = M(i,:);
    M(1,:) = M(1,:)*sign(R(1,1)); M(2,:) = M(2,:)*sign(R(2,2));
    PFMts{s} = M; clear M
end

% Add PROFUMO results
PFMmaps_matrix = zeros(V,D,S);
PFMcorrs = zeros(s,3);
for s = 1:S
    d = Y((s-1)*t+1:s*t,:)';
    NODEts = PFMts{s}';
    MAP = PFMmaps{s};
    PFMmaps_matrix(:,:,s) = MAP;
    PFMcorrs(s,1) = corr(NODEts(:,1),NODEts(:,2)); 
    R = PFMnet{s}; PFMcorrs(s,2) = R(2,1); clear R
    PFMcorrs(s,3) = corr(MAP(:,1),MAP(:,2)); 
    MapTcorrs(s,5) = corr(MAP(:,1),GTmaps(:,(s-1)*D+1));
    MapTcorrs(s,6) = corr(MAP(:,2),GTmaps(:,(s-1)*D+2));
    MapTcorrs(s,7) = corr(NODEts(:,1),GTts(:,(s-1)*D+1));
    MapTcorrs(s,8) = corr(NODEts(:,2),GTts(:,(s-1)*D+2));
end

% Plot Results
figure; set(gcf,'Position',[1150 450 800 600],'PaperPositionMode','auto')
M = [GTcorrs DRcorrs PFMcorrs MapTcorrs];
imagesc(corr(M),[-1 1]); colorbar
hline([2.5 4.5 7.5],'k'); vline([2.5 4.5 7.5],'k');
N = {'GT Tnets','GT Snets','DR Tnets','DR Snets',...
    'PFM Tnets1','PFM Tnets2','PFM Snets',...
    'Map1 corr','Map2 corr','TS1 corr','TS2 corr'};
set(gca,'xtick',1:length(N),'xticklabel',N,'ytick',1:length(N),'yticklabel',N,'FontSize',14)
for n = 1:size(M,2); text(7.7,n,sprintf('%1.4f',mean(M(:,n))),'FontSize',14); end
title('correlations across subjects','FontSize',16);

figure; set(gcf,'Position',[150 75 900 1000],'PaperPositionMode','auto')
M = [mean(reshape(GTmaps,V,D,S),3) icaS' mean(DRmaps,3) PFMgmap];
M = M.*repmat(sign(mean(M)),size(M,1),1)./repmat(max(abs(M)),size(M,1),1);
imagesc(M);
vline([2.5 4.5],'k'); colorbar
R = corr(mean(reshape(GTmaps,V,D,S),3)); text(1,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
R = corr(icaS'); text(3,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
R = corr(mean(DRmaps,3)); text(5,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
R = corr(PFMgmap); text(7,9500,sprintf('r=%1.2f',R(2,1)),'FontSize',14);
N = {'GT map1','GT map2','gICA map1','gICA map2','DR map1','DR map2','PFM map1','PFM map2'};
set(gca,'xtick',1:length(M),'xticklabel',N,'FontSize',14)
title('Group average maps','FontSize',16)
