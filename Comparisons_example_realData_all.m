clear all; close all; clc

% Find matching components
cifti = ft_read_cifti('../PFM/maps_full_originalOrder.dtseries.nii'); gPFM = cifti.dtseries; dtseries_isnan = isnan(gPFM(:,1)); gPFM(isnan(gPFM(:,1))==1,:) = [];
gICA = ft_read_cifti(fullfile('/vols/Data/HCP/Phase2/group1200/groupICA/','groupICA_3T_HCP1200_MSMAll_d50.ica','melodic_IC.dtseries.nii')); gICA = gICA.dtseries; gICA(isnan(gICA(:,1))==1,:) = [];
[C12,munkres_assign] = spatialcorr(gICA,gPFM);
for n = 1:size(C12,1); J50(n) = find(munkres_assign(n,:)==1); M(n) = C12(n,J50(n)); end
[~,I50] = sort(abs(M),'descend');
J50 = J50(I50);
J50 = [J50 setdiff(1:size(munkres_assign,2),J50)];
SpatCor50 = C12(I50,J50); SpatCor50 = SpatCor50(eye(50)==1); L = find(SpatCor50>0.7,1,'last');
Iica = I50(1:L); Ipfm = J50(1:L);

addpath /vols/Scratch/janineb/HCP/DMN/DMN_functions/hline_vline/
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
HCPdir = '/vols/Data/HCP/Phase2/subjects1200';
INdir = '/vols/Scratch/HCP/rfMRI/PROFUMO/S1003_M50.pfm/FinalModel';
scans = {'rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR','rfMRI_REST2_RL'};
Pdir = '/vols/Scratch/HCP/rfMRI/PROFUMO/S1003_M50.pfm/FinalModel/';
subs = dir('/vols/Scratch/HCP/rfMRI/PROFUMO/S820_M50_Aug16.pfm/FinalModel/Subjects/');
subs = subs(3:end); subs(169) = [];

% Group maps:
gICA = gICA(:,Iica); 
gICA_thres = gICA; gICA_thres(abs(gICA_thres)<prctile(abs(gICA(:)),95)) = 0;
gPFM = gPFM(:,Ipfm); 

PFMall = [];
ICAall = [];
I = 1;
for x = 1:length(Iica);
    for y = 1:length(Iica);
        if x>y
            if abs(corr(gPFM(:,x),gPFM(:,y)))>0.3
                PFMall = [PFMall gPFM(:,[x y])]; ICAall = [ICAall gICA(:,[x y])];
                [Ms1,Is1] = sort(mean([gPFM(:,x) gICA(:,x)],2),'descend');
                [Ms2,Is2] = sort(mean([gPFM(:,y) gICA(:,y)],2),'descend');
                %[Ms1,Is1] = sort(gPFM(:,x),'descend');
                %[Ms2,Is2] = sort(gPFM(:,y),'descend');
                Is = [Is1(1:find(Ms1<0.2,1,'first'))' Is2(1:find(Ms2<0.2,1,'first'))']; Is = [Is setdiff(1:91282,Is)];
                Grays = 1:find(Ms1<0.2,1,'first')+find(Ms2<0.2,1,'first');
                figure; imagesc([gPFM(Is(Grays),[x y]) gICA(Is(Grays),[x y])],[-1 1]);
                set(gcf,'Position',[0 0 900 945],'PaperPositionMode','auto')
                vline(2.5,'k')
                set(gca,'xtick',1:4,'xticklabel',{sprintf('PFM %d',Ipfm(x)) sprintf('PFM %d',Ipfm(y)) sprintf('ICA %d',Iica(x)) sprintf('ICA %d',Iica(y))})
                title(sprintf('PFMr = %1.2f, ICAr = %1.2f, ICAr_thres = %1.2f',corr(gPFM(:,x),gPFM(:,y)),corr(gICA(:,x),gICA(:,y)),corr(gICA_thres(:,x),gICA_thres(:,y))));
                print(sprintf('%02d_example_real',I),'-dpng','-r150'); I = I+1;
            end
        end
    end
end
%save('temp_real.mat','ICAall','PFMall')
