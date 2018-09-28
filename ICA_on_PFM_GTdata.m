clear all; close all; clc
% Indices (matched maps based on Spatial_correlation_ICA.m):
Iica = [1 23 2 20 4 3 5 8 11 19 10];
Ipfm = [1 2 12 3 7 11 22 10 25 15 18];

addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
addpath /vols/Scratch/janineb/HCP/netmat_simulations/
addpath /vols/Scratch/janineb/matlab/FastICA_25/
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
HCPdir = '/vols/Data/HCP/Phase2/subjects1200';
INdir = '/vols/Scratch/HCP/rfMRI/PROFUMO/S1003_M50.pfm/FinalModel';
scans = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
Pdir = '/vols/Scratch/HCP/rfMRI/PROFUMO/S1003_M50.pfm/FinalModel/';
subs = dir('/vols/Scratch/HCP/rfMRI/PROFUMO/S820_M50_Aug16.pfm/FinalModel/Subjects/');
subs = subs(3:end); subs(169) = [];

% Group maps:
groupICA = '/vols/Data/HCP/Phase2/group1200/groupICA/groupICA_3T_HCP1200_MSMAll_d25.ica/melodic_IC.dtseries.nii';
gICA = ft_read_cifti(groupICA); gICA = gICA.dtseries(:,Iica); gICA(isnan(gICA(:,1))==1,:) = [];
gICA_thres = gICA; gICA_thres(abs(gICA)<0.1) = 0;
groupPFM = '../PFM/maps_full_originalOrder_max.dtseries.nii';
gPFM = ft_read_cifti(groupPFM); gPFM = gPFM.dtseries(:,Ipfm); gPFM(isnan(gPFM(:,1))==1,:) = [];
%gICA_norm = gICA.*repmat(sign(mean(gPFM)),size(gICA,1),1)./repmat(max(abs(gICA)),size(gICA,1),1);
%gPFM_norm = gPFM.*repmat(sign(mean(gPFM)),size(gPFM,1),1)./repmat(max(abs(gPFM)),size(gPFM,1),1);

% Subject maps:
s = 2; runs = {'1_LR','1_RL','2_LR','2_RL'};
means = h5read(fullfile(Pdir,'Subjects',subs(s).name,'SpatialMaps.post','Signal','Means.hdf5'),'/dataset');
probs = h5read(fullfile(Pdir,'Subjects',subs(s).name,'SpatialMaps.post','Signal','MembershipProbabilities.hdf5'),'/dataset');
sPFM = means .* probs; sPFM = sPFM(:,Ipfm); clear means probs
TC = h5read(fullfile(Pdir,'Subjects',subs(s).name,'Runs',runs{1},'TimeCourses.post','CleanTimeCourses','Means.hdf5'),'/dataset')';
a = importfile(fullfile(Pdir,'Subjects', subs(s).name,'Runs',runs{1},'NoisePrecision.post','GammaPosterior.txt'));
noiseStd = sqrt(a(2)/a(1));
H = h5read(fullfile(Pdir,'Subjects',subs(s).name,'Runs',runs{1},'ComponentWeightings.post','Means.hdf5'),'/dataset');
TC = TC(:,Ipfm); H = H(Ipfm);

PFMall = [];
ICAall = [];
I = 1;
for x = 1:length(Iica);
    for y = 1:length(Iica);
        if x>y
            if abs(corr(gPFM(:,x),gPFM(:,y)))>0.2
                if abs(corr(sPFM(:,x),sPFM(:,y)))>0.2
                    
                    D = sPFM(:,[x y]) * diag(H([x y])) * TC(:,[x y])';
                    [sICA,A,W] = fastica(D','approach','symm','g','tanh','lastEig',2);
                    sICA = sICA';
                    R = corr([sPFM(:,[x y]) sICA]); R = R(3:end,1:2);
                    [~,i] = max(abs(R)); if length(unique(i))==1; noi = setdiff(1:2,i); [~,inew] = max(abs(R(noi,:))); i(inew) = noi; end
                    sICA(:,1) = sICA(:,1)*sign(R(1,i(1))); A(:,1) = A(:,1)*sign(R(1,i(1)));
                    sICA(:,2) = sICA(:,2)*sign(R(2,i(2))); A(:,2) = A(:,2)*sign(R(2,i(2)));
                    sICA = sICA(:,i);
                    sICA_thres = sICA; sICA_thres(abs(sICA_thres)<prctile(abs(sICA(:)),95)) = 0;
                    PFMall = [PFMall sPFM(:,[x y])]; ICAall = [ICAall sICA];
                    
                    [Ms1,Is1] = sort(mean([sPFM(:,x) sICA(:,1)],2),'descend');
                    [Ms2,Is2] = sort(mean([sPFM(:,y) sICA(:,2)],2),'descend');
                    Is = [Is1(1:find(Ms1<0.2,1,'first'))' Is2(1:find(Ms2<0.2,1,'first'))']; Is = [Is setdiff(1:91282,Is)];
                    Grays = 1:find(Ms1<0.2,1,'first')+find(Ms2<0.2,1,'first');
                    figure; imagesc([15*gPFM(Is(Grays),[x y]) sPFM(Is(Grays),[x y]) sICA(Is(Grays),:)],[-15 15]);
                    text(1,1000,'Group PFM','FontSize',16);
                    text(3,1000,'Subject PFM','FontSize',16);
                    text(5,1000,'Subject ICA','FontSize',16);
                    set(gcf,'Position',[0 0 900 945],'PaperPositionMode','auto')
                    vline([2.5 4.5],'k')
                    set(gca,'xtick',2:6,'xticklabel',{sprintf('%1.2f | %1.2f',corr(TC(:,x),TC(:,y)),corr(A(:,1),A(:,2))) sprintf('%1.2f',corr(gPFM(:,x),sPFM(:,x))) sprintf('%1.2f',corr(gPFM(:,y),sPFM(:,y))) sprintf('%1.2f | %1.2f',corr(gPFM(:,x),sICA(:,1)),corr(sPFM(:,x),sICA(:,1))) sprintf('%1.2f | %1.2f',corr(gPFM(:,y),sICA(:,2)),corr(sPFM(:,y),sICA(:,2)))})
                    title(sprintf('PFMg = %1.2f, PFMs = %1.2f, ICAs = %1.2f, ICAs thres = %1.2f',corr(gPFM(:,x),gPFM(:,y)),corr(sPFM(:,x),sPFM(:,y)),corr(sICA(:,1),sICA(:,2)),corr(sICA_thres(:,1),sICA_thres(:,2))));
                    print(sprintf('%02d_example',I),'-dpng','-r150'); I = I+1;
                end
            end
        end
    end
end
save('temp.mat','ICAall','PFMall')
