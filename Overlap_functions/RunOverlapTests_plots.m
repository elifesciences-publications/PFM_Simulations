

pfmPg = pfmPg1;
[C12,munkres_assign] = spatialcorr(pfmPg,Pg);
[i_pfm,j] = find(munkres_assign==1); sign_pfm = sign(C12(munkres_assign==1));
pfmPg_new = pfmPg1_new;
[C12,munkres_assign] = spatialcorr(pfmPg_new,Pg);
[i_pfm_new,j] = find(munkres_assign==1); sign_pfm_new = sign(C12(munkres_assign==1));
sicaPg = sicaP1{1};
[C12,munkres_assign] = spatialcorr(sicaPg,Pg);
[i_ica,j] = find(munkres_assign==1); sign_ica = sign(C12(munkres_assign==1));
[C12,munkres_assign] = spatialcorr(sicaPg_new,Pg);
[i_ica_new,j] = find(munkres_assign==1); sign_ica_new = sign(C12(munkres_assign==1));
figure; set(gcf,'Position',[900 570 1000 415],'PaperPositionMode','auto')
subplot(1,5,1); imagesc(Pg); title('ground truth group maps')
subplot(1,5,2); imagesc(pfmPg(:,i_pfm).*repmat(sign_pfm',atlasParams.V,1)); title('PFM old group maps'); colormap parula
subplot(1,5,3); imagesc(pfmPg_new(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1)); title('PFM new group maps'); colormap parula
subplot(1,5,4); imagesc(sicaPg(:,i_ica).*repmat(sign_ica',atlasParams.V,1)); title('sICA group maps'); colormap parula
subplot(1,5,5); imagesc(sicaPg_new(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1)); title('sICA melodic group maps'); colormap parula
GT_Snet = zeros(params.iN,params.iN,params.S);
GT_Tnet = zeros(params.iN,params.iN,params.S);
PFM_Snet = zeros(params.iN,params.iN,params.S);
PFM_Tnet = zeros(params.iN,params.iN,params.S);
PFMnew_Snet = zeros(params.iN,params.iN,params.S);
PFMnew_Tnet = zeros(params.iN,params.iN,params.S);
sICA_Snet = zeros(params.iN,params.iN,params.S);
sICA_Tnet = zeros(params.iN,params.iN,params.S);
sICAnew_Snet = zeros(params.iN,params.iN,params.S);
sICAnew_Tnet = zeros(params.iN,params.iN,params.S);

for s = 1:params.S
    GT_Snet(:,:,s) = corr(P{s});
    GT_Tnet(:,:,s) = corr(Ag{s}{1}');
    p = pfmP1{s}; PFM_Snet(:,:,s) = corr(p(:,i_pfm).*repmat(sign_pfm',atlasParams.V,1)); clear p
    a = pfmA1{s}{1}'; PFM_Tnet(:,:,s) = corr(a(:,i_pfm).*repmat(sign_pfm',params.T,1)); clear a
    p = pfmP1_new{s}; PFMnew_Snet(:,:,s) = corr(p(:,i_pfm_new).*repmat(sign_pfm_new',atlasParams.V,1)); clear p
    a = pfmA1_new{s}{1}'; PFMnew_Tnet(:,:,s) = corr(a(:,i_pfm_new).*repmat(sign_pfm_new',params.T,1)); clear a
    p = sica_P1_DRtest{s}; sICA_Snet(:,:,s) = corr(p(:,i_ica).*repmat(sign_ica',atlasParams.V,1)); clear p
    a = sica_A1_DRtest{s}{1}; sICA_Tnet(:,:,s) = corr(a(:,i_ica).*repmat(sign_ica',params.T,1)); clear a
    p = sica_P1_DR_new{s}; sICAnew_Snet(:,:,s) = corr(p(:,i_ica_new).*repmat(sign_ica_new',atlasParams.V,1)); clear p
    a = sica_A1_DR_new{s}{1}; sICAnew_Tnet(:,:,s) = corr(a(:,i_ica_new).*repmat(sign_ica_new',params.T,1)); clear a
end

% Plot subject simlarity of netmats
A = ones(params.iN); A = triu(A,1); A = find(A);
Pt = reshape(PFM_Tnet,params.iN*params.iN,params.S); Pt_new = reshape(PFMnew_Tnet,params.iN*params.iN,params.S);
It = reshape(sICA_Tnet,params.iN*params.iN,params.S); It_new = reshape(sICAnew_Tnet,params.iN*params.iN,params.S);
Gt = reshape(GT_Tnet,params.iN*params.iN,params.S);
figure; set(gcf,'Position',[900 0 1000 455],'PaperPositionMode','auto')
R = corr([Gt(A,:) Pt(A,:) Pt_new(A,:) It(A,:) It_new(A,:)]);
imagesc(R,[-1 1]); colorbar; colormap parula
hline([30.5 60.5 90.5 120.5],'k'); vline([30.5 60.5 90.5 120.5],'k')  
text(5,params.S+2,'ground truth - PFM old'); r = R(1:params.S,params.S+1:params.S*2); r = r(eye(params.S)==1); text(5,params.S*2-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*2+2,'ground truth - PFM new'); r = R(1:params.S,params.S*2+1:params.S*3); r = r(eye(params.S)==1); text(5,params.S*3-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*3+2,'ground truth - sICA'); r = R(1:params.S,params.S*3+1:params.S*4); r = r(eye(params.S)==1); text(5,params.S*4-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*4+2,'ground truth - sICA melodic'); r = R(1:params.S,params.S*4+1:params.S*5); r = r(eye(params.S)==1); text(5,params.S*5-5,sprintf('mean r=%1.2f',mean(r)));
set(gca,'ytick',15:30:params.S*5,'yticklabel',{'ground truth','PFM old','PFM new','sICA','sICA melodic'})
title('Subject-level similarities of full netmats')
  
% Plot group means
figure; set(gcf,'Position',[0 0 900 955],'PaperPositionMode','auto')
subplot(5,2,1); imagesc(mean(GT_Snet,3),[-0.5 0.5]); colorbar; title('Ground truth spatial correlations (group mean)'); colormap parula
subplot(5,2,2); imagesc(mean(GT_Tnet,3),[-0.5 0.5]); colorbar; title('Ground truth temporal correlations (group mean)'); colormap parula
subplot(5,2,3); imagesc(mean(PFM_Snet,3),[-0.5 0.5]); colorbar; title('PFM old spatial correlations (group mean)'); colormap parula
subplot(5,2,4); imagesc(mean(PFM_Tnet,3),[-0.5 0.5]); colorbar; title('PFM old temporal correlations (group mean)'); colormap parula
subplot(5,2,5); imagesc(mean(PFMnew_Snet,3),[-0.5 0.5]); colorbar; title('PFM new spatial correlations (group mean)'); colormap parula
subplot(5,2,6); imagesc(mean(PFMnew_Tnet,3),[-0.5 0.5]); colorbar; title('PFM new temporal correlations (group mean)'); colormap parula
subplot(5,2,7); imagesc(mean(sICA_Snet,3),[-0.5 0.5]); colorbar; title('sICA DR spatial correlations (group mean)'); colormap parula
subplot(5,2,8); imagesc(mean(sICA_Tnet,3),[-0.5 0.5]); colorbar; title('sICA DR temporal correlations (group mean)'); colormap parula
subplot(5,2,9); imagesc(mean(sICAnew_Snet,3),[-0.5 0.5]); colorbar; title('sICA melodic DR spatial correlations (group mean)'); colormap parula
subplot(5,2,10); imagesc(mean(sICAnew_Tnet,3),[-0.5 0.5]); colorbar; title('sICA melodic DR temporal correlations (group mean)'); colormap parula

% Plot scatter plots
AllNets = cat(3,reshape(GT_Tnet,params.iN*params.iN,params.S)',reshape(PFM_Tnet,params.iN*params.iN,params.S)',reshape(sICAnew_Tnet,params.iN*params.iN,params.S)',...
    reshape(GT_Snet,params.iN*params.iN,params.S)',reshape(PFM_Snet,params.iN*params.iN,params.S)',reshape(sICAnew_Snet,params.iN*params.iN,params.S)');
AllNames = {'GT temporal net','PFM old temporal net','sICA melodic temporal net','GT spatial net','PFM old spatial net','sICA melodic spatial net'};
A = ones(params.iN); A = triu(A,1); A = find(A==1); An = length(A);
figure; set(gcf,'Position',[0 0 1900 955],'PaperPositionMode','auto')
ns = [1 1 1 1 4 4 4 4 2 3]; ms = [2 3 5 6 5 6 2 3 5 6];
for i = 1:length(ns);
    subplot(3,4,i);
    scatter(reshape(AllNets(:,A,ns(i)),params.S*An,1),reshape(AllNets(:,A,ms(i)),params.S*An,1)); axis([-1 1 -1 1])
    xlabel(AllNames{ns(i)}); ylabel(AllNames{ms(i)});
    title(sprintf('r = %1.2f',corr(reshape(AllNets(:,A,ns(i)),params.S*An,1),reshape(AllNets(:,A,ms(i)),params.S*An,1))));
end

%%% SAME PLOTS BUT REMOVING EMPTY PFM maps
PFMempty = find(sum(abs(pfmPg_new))<10); PFMgood = setdiff(1:params.iN,PFMempty); Ngood = length(PFMgood);

% Plot subject simlarity of netmats
A = ones(Ngood); A = triu(A,1); A = find(A);
Pt = reshape(PFM_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S); Pt_new = reshape(PFMnew_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S);
It = reshape(sICA_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S); It_new = reshape(sICAnew_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S);
Gt = reshape(GT_Tnet(PFMgood,PFMgood,:),Ngood*Ngood,params.S);
figure; set(gcf,'Position',[900 0 1000 455],'PaperPositionMode','auto')
R = corr([Gt(A,:) Pt(A,:) Pt_new(A,:) It(A,:) It_new(A,:)]);
imagesc(R,[-1 1]); colorbar; colormap parula
hline([30.5 60.5 90.5 120.5],'k'); vline([30.5 60.5 90.5 120.5],'k')  
text(5,params.S+2,'ground truth - PFM old'); r = R(1:params.S,params.S+1:params.S*2); r = r(eye(params.S)==1); text(5,params.S*2-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*2+2,'ground truth - PFM new'); r = R(1:params.S,params.S*2+1:params.S*3); r = r(eye(params.S)==1); text(5,params.S*3-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*3+2,'ground truth - sICA'); r = R(1:params.S,params.S*3+1:params.S*4); r = r(eye(params.S)==1); text(5,params.S*4-5,sprintf('mean r=%1.2f',mean(r)));
text(5,params.S*4+2,'ground truth - sICA melodic'); r = R(1:params.S,params.S*4+1:params.S*5); r = r(eye(params.S)==1); text(5,params.S*5-5,sprintf('mean r=%1.2f',mean(r)));
set(gca,'ytick',15:30:params.S*5,'yticklabel',{'ground truth','PFM old','PFM new','sICA','sICA melodic'})
title('Subject-level similarities of full netmats')
  
% Plot group means
figure; set(gcf,'Position',[0 0 900 955],'PaperPositionMode','auto')
subplot(5,2,1); imagesc(mean(GT_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('Ground truth spatial correlations (group mean)'); colormap parula
subplot(5,2,2); imagesc(mean(GT_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('Ground truth temporal correlations (group mean)'); colormap parula
subplot(5,2,3); imagesc(mean(PFM_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM old spatial correlations (group mean)'); colormap parula
subplot(5,2,4); imagesc(mean(PFM_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM old temporal correlations (group mean)'); colormap parula
subplot(5,2,5); imagesc(mean(PFMnew_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM new spatial correlations (group mean)'); colormap parula
subplot(5,2,6); imagesc(mean(PFMnew_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('PFM new temporal correlations (group mean)'); colormap parula
subplot(5,2,7); imagesc(mean(sICA_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA DR spatial correlations (group mean)'); colormap parula
subplot(5,2,8); imagesc(mean(sICA_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA DR temporal correlations (group mean)'); colormap parula
subplot(5,2,9); imagesc(mean(sICAnew_Snet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA melodic DR spatial correlations (group mean)'); colormap parula
subplot(5,2,10); imagesc(mean(sICAnew_Tnet(PFMgood,PFMgood,:),3),[-0.5 0.5]); colorbar; title('sICA melodic DR temporal correlations (group mean)'); colormap parula

