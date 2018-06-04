addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/matlab/

pfmPg = pfmPg1;
[C12,munkres_assign] = spatialcorr(pfmPg,Pg);
[i_pfm,j] = find(munkres_assign==1); sign_pfm = sign(C12(munkres_assign==1));
pfmPg_new = pfmPg1_new;
[C12,munkres_assign] = spatialcorr(pfmPg_new,Pg);
[i_pfm_new,j] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
ticaPg = ticaP1{1};
[C12,munkres_assign] = spatialcorr(ticaPg,Pg);
[i_tica,j] = find(munkres_assign==1); sign_tica = sign(C12(munkres_assign==1));
[C12,munkres_assign] = spatialcorr(sicaPg_new,Pg);
[i_ica_new,j] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
figure; set(gcf,'Position',[900 570 1000 415],'PaperPositionMode','auto')
subplot(1,5,1); imagesc(Pg); title('ground truth group maps')
subplot(1,5,2); imagesc(pfmPg(:,i_pfm).*repmat(sign_pfm',atlasParams.V,1)); title('PFM old group maps'); colormap parula
subplot(1,5,3); imagesc(pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1)); title('PFM new group maps'); colormap parula
subplot(1,5,4); imagesc(ticaPg(:,i_tica).*repmat(sign_tica',atlasParams.V,1)); title('tICA group maps'); colormap parula
subplot(1,5,5); imagesc(sicaPg_new(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1)); title('sICA melodic group maps'); colormap parula
GT_Snet = zeros(params.iN,params.iN,params.S);
GT_Tnet = zeros(params.iN,params.iN,params.S);
PFM_Snet = zeros(params.iN,params.iN,params.S);
PFM_Tnet = zeros(params.iN,params.iN,params.S);
PFMnew_Snet = zeros(params.iN,params.iN,params.S);
PFMnew_Tnet = zeros(params.iN,params.iN,params.S);
tICA_Snet = zeros(params.iN,params.iN,params.S);
tICA_Tnet = zeros(params.iN,params.iN,params.S);
sICAnew_Snet = zeros(params.iN,params.iN,params.S);
sICAnew_Tnet = zeros(params.iN,params.iN,params.S);

for s = 1:params.S
    GT_Snet(:,:,s) = corr(P{s});
    GT_Tnet(:,:,s) = corr(Ag{s}{1}');
    p = pfmP1{s}; PFM_Snet(:,:,s) = corr(p(:,i_pfm).*repmat(sign_pfm',atlasParams.V,1)); clear p
    a = pfmA1{s}{1}'; PFM_Tnet(:,:,s) = corr(a(:,i_pfm).*repmat(sign_pfm',params.T,1)); clear a
    p = pfmP1_new{s}; PFMnew_Snet(:,:,s) = corr(p(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1)); clear p
    a = pfmA1_new{s}{1}'; PFMnew_Tnet(:,:,s) = corr(a(:,i_pfm_new).*repmat(sign_pfm_new',params.T,1)); clear a
    p = ticaP1_DR{s}; tICA_Snet(:,:,s) = corr(p(:,i_tica).*repmat(sign_tica',atlasParams.V,1)); clear p
    a = ticaA1_DR{s}{1}'; tICA_Tnet(:,:,s) = corr(a(:,i_tica).*repmat(sign_tica',params.T,1)); clear a
    p = sica_P1_DR_new{s}; sICAnew_Snet(:,:,s) = corr(p(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1)); clear p
    a = sica_A1_DR_new{s}{1}; sICAnew_Tnet(:,:,s) = corr(a(:,i_ica_new).*repmat(sign_ica_new',params.T,1)); clear a
end

% Plot subject simlarity of netmats
A = ones(params.iN); A = triu(A,1); A = find(A);
Pt = reshape(PFM_Tnet,params.iN*params.iN,params.S); Pt_new = reshape(PFMnew_Tnet,params.iN*params.iN,params.S);
It = reshape(tICA_Tnet,params.iN*params.iN,params.S); It_new = reshape(sICAnew_Tnet,params.iN*params.iN,params.S);
Gt = reshape(GT_Tnet,params.iN*params.iN,params.S);
figure; set(gcf,'Position',[900 0 1000 455],'PaperPositionMode','auto')
R = corr([Gt(A,:) Pt(A,:) Pt_new(A,:) It(A,:) It_new(A,:)]);
imagesc(R,[-1 1]); colorbar; colormap parula
hline([30.5 60.5 90.5 120.5],'k'); vline([30.5 60.5 90.5 120.5],'k')  
text(5,params.S+2,'ground truth - PFM old'); r = R(1:params.S,params.S+1:params.S*2); r = r(eye(params.S)==1); text(5,params.S*2-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*2+2,'ground truth - PFM new'); r = R(1:params.S,params.S*2+1:params.S*3); r = r(eye(params.S)==1); text(5,params.S*3-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*3+2,'ground truth - tICA'); r = R(1:params.S,params.S*3+1:params.S*4); r = r(eye(params.S)==1); text(5,params.S*4-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*4+2,'ground truth - sICA melodic'); r = R(1:params.S,params.S*4+1:params.S*5); r = r(eye(params.S)==1); text(5,params.S*5-5,sprintf('mean r=%1.2f',mean(r)));
set(gca,'ytick',15:30:params.S*5,'yticklabel',{'ground truth','PFM old','PFM new','tICA','sICA melodic'})
title('Subject-level similarities of full netmats')
  
% Plot group means
figure; set(gcf,'Position',[0 0 900 955],'PaperPositionMode','auto')
subplot(5,2,1); imagesc(mean(GT_Snet,3),[-0.5 0.5]); colorbar; title('Ground truth spatial correlations (group mean)'); colormap parula
subplot(5,2,2); imagesc(mean(GT_Tnet,3),[-0.5 0.5]); colorbar; title('Ground truth temporal correlations (group mean)'); colormap parula
subplot(5,2,3); imagesc(mean(PFM_Snet,3),[-0.5 0.5]); colorbar; title('PFM old spatial correlations (group mean)'); colormap parula
subplot(5,2,4); imagesc(mean(PFM_Tnet,3),[-0.5 0.5]); colorbar; title('PFM old temporal correlations (group mean)'); colormap parula
subplot(5,2,5); imagesc(mean(PFMnew_Snet,3),[-0.5 0.5]); colorbar; title('PFM new spatial correlations (group mean)'); colormap parula
subplot(5,2,6); imagesc(mean(PFMnew_Tnet,3),[-0.5 0.5]); colorbar; title('PFM new temporal correlations (group mean)'); colormap parula
subplot(5,2,7); imagesc(mean(tICA_Snet,3),[-0.5 0.5]); colorbar; title('tICA DR spatial correlations (group mean)'); colormap parula
subplot(5,2,8); imagesc(mean(tICA_Tnet,3),[-0.5 0.5]); colorbar; title('tICA DR temporal correlations (group mean)'); colormap parula
subplot(5,2,9); imagesc(mean(sICAnew_Snet,3),[-0.5 0.5]); colorbar; title('sICA melodic DR spatial correlations (group mean)'); colormap parula
subplot(5,2,10); imagesc(mean(sICAnew_Tnet,3),[-0.5 0.5]); colorbar; title('sICA melodic DR temporal correlations (group mean)'); colormap parula

% Plot scatter plots
AllNets = cat(3,reshape(GT_Tnet,params.iN*params.iN,params.S)',reshape(PFMnew_Tnet,params.iN*params.iN,params.S)',reshape(sICAnew_Tnet,params.iN*params.iN,params.S)',...
    reshape(GT_Snet,params.iN*params.iN,params.S)',reshape(PFMnew_Snet,params.iN*params.iN,params.S)',reshape(sICAnew_Snet,params.iN*params.iN,params.S)');
AllNames = {'GT temporal net','PFM new temporal net','sICA melodic temporal net','GT spatial net','PFM new spatial net','sICA melodic spatial net'};
A = ones(params.iN); A = triu(A,1); A = find(A==1); An = length(A);
figure; set(gcf,'Position',[0 0 1900 955],'PaperPositionMode','auto')
ns = [1 1 1 1 4 4 4 4 2 3]; ms = [2 3 5 6 5 6 2 3 5 6];
for i = 1:length(ns);
    subplot(3,4,i);
    scatplot(reshape(AllNets(:,A,ns(i)),params.S*An,1),reshape(AllNets(:,A,ms(i)),params.S*An,1)); 
    axis([-1 1 -1 1])
    xlabel(AllNames{ns(i)}); ylabel(AllNames{ms(i)});
    r = corr(reshape(AllNets(:,A,ns(i)),params.S*An,1),reshape(AllNets(:,A,ms(i)),params.S*An,1));
    title(sprintf('r = %1.2f',r));
    hold on; if r>0; plot(-1:0.1:1,-1:0.1:1,'r'); else; plot(-1:0.1:1,1:-0.1:-1,'r'); end
end

%%% PLOT INDIVIDUAL MAPS
Pg_norm = Pg;
Pg_norm = Pg_norm.*repmat(sign(mean(Pg_norm)),size(Pg_norm,1),1)./repmat(max(abs(Pg_norm)),size(Pg_norm,1),1);
pfmPg_new_norm = pfmPg_new;
pfmPg_new_norm = pfmPg_new_norm.*repmat(sign(mean(pfmPg_new_norm)),size(pfmPg_new_norm,1),1)./repmat(max(abs(pfmPg_new_norm)),size(pfmPg_new_norm,1),1);
sicaPg_new_norm = sicaPg_new;
sicaPg_new_norm = sicaPg_new_norm.*repmat(sign(mean(sicaPg_new_norm)),size(sicaPg_new_norm,1),1)./repmat(max(abs(sicaPg_new_norm)),size(sicaPg_new_norm,1),1);
ticaPg_new_norm = ticaPg;
ticaPg_new_norm = ticaPg_new_norm.*repmat(sign(mean(ticaPg_new_norm)),size(ticaPg_new_norm,1),1)./repmat(max(abs(ticaPg_new_norm)),size(ticaPg_new_norm,1),1);
figure; set(gcf,'Position',[0 0 1900 955],'PaperPositionMode','auto')
[~,M] = max(sum(abs(Pg),2));
for n = 1:15
    subplot(3,5,n);
    plot(sicaPg_new_norm(:,i_ica_new(n)),'g','linewidth',2); 
    hold on
    plot(ticaPg_new_norm(:,i_ica_new(n)),'m','linewidth',2);
    plot(pfmPg_new_norm(:,i_pfm_new(n)),'r','linewidth',2);
    plot(Pg_norm(:,n),'linewidth',2); 
    vline(M); hold off
end
legend({'sICA new','tICA','PFM new','ground truth'})

% %%% SAME PLOTS BUT REMOVING EMPTY PFM maps
% PFMempty = find(sum(abs(pfmPg_new))<10); PFMgood = setdiff(1:params.iN,PFMempty); Ngood = length(PFMgood);
% 
% % Plot subject simlarity of netmats
% A = ones(Ngood); A = triu(A,1); A = find(A);
% Pt = reshape(PFM_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S); Pt_new = reshape(PFMnew_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S);
% It = reshape(sICA_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S); It_new = reshape(sICAnew_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S);
% Gt = reshape(GT_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S);
% figure; set(gcf,'Position',[900 0 1000 455],'PaperPositionMode','auto')
% R = corr([Gt(A,:) Pt(A,:) Pt_new(A,:) It(A,:) It_new(A,:)]);
% imagesc(R,[-1 1]); colorbar; colormap parula
% hline([30.5 60.5 90.5 120.5],'k'); vline([30.5 60.5 90.5 120.5],'k')  
% text(5,params.S+2,'ground truth - PFM old'); r = R(1:params.S,params.S+1:params.S*2); r = r(eye(params.S)==1); text(5,params.S*2-5,sprintf('mean r=%1.2f',mean(r)));
% text(5,params.S*2+2,'ground truth - PFM new'); r = R(1:params.S,params.S*2+1:params.S*3); r = r(eye(params.S)==1); text(5,params.S*3-5,sprintf('mean r=%1.2f',mean(r)));
% text(5,params.S*3+2,'ground truth - sICA'); r = R(1:params.S,params.S*3+1:params.S*4); r = r(eye(params.S)==1); text(5,params.S*4-5,sprintf('mean r=%1.2f',mean(r)));
% text(5,params.S*4+2,'ground truth - sICA melodic'); r = R(1:params.S,params.S*4+1:params.S*5); r = r(eye(params.S)==1); text(5,params.S*5-5,sprintf('mean r=%1.2f',mean(r)));
% set(gca,'ytick',15:30:params.S*5,'yticklabel',{'ground truth','PFM old','PFM new','sICA','sICA melodic'})
% title('Subject-level similarities of full netmats')
%   
% % Plot group means
% figure; set(gcf,'Position',[0 0 900 955],'PaperPositionMode','auto')
% subplot(5,2,1); imagesc(mean(GT_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('Ground truth spatial correlations (group mean)'); colormap parula
% subplot(5,2,2); imagesc(mean(GT_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('Ground truth temporal correlations (group mean)'); colormap parula
% subplot(5,2,3); imagesc(mean(PFM_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM old spatial correlations (group mean)'); colormap parula
% subplot(5,2,4); imagesc(mean(PFM_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM old temporal correlations (group mean)'); colormap parula
% subplot(5,2,5); imagesc(mean(PFMnew_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM new spatial correlations (group mean)'); colormap parula
% subplot(5,2,6); imagesc(mean(PFMnew_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM new temporal correlations (group mean)'); colormap parula
% subplot(5,2,7); imagesc(mean(sICA_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA DR spatial correlations (group mean)'); colormap parula
% subplot(5,2,8); imagesc(mean(sICA_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA DR temporal correlations (group mean)'); colormap parula
% subplot(5,2,9); imagesc(mean(sICAnew_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA melodic DR spatial correlations (group mean)'); colormap parula
% subplot(5,2,10); imagesc(mean(sICAnew_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA melodic DR temporal correlations (group mean)'); colormap parula
