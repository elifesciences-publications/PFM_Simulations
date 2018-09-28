clear all; close all; clc

% Set paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/saveJSONfile/')
path(path,'/vols/Scratch/janineb/HCP/DMN1200/Functions')

% Load parameters
load('Results/ICA_DRtest.mat');
%load('Results/ICA_DRtest_NoSpatCorr.mat');
%load('Results/ICA_DRtest_NoTempCorr.mat');

% Save results as niftis ready for PROFUMO
TR = 0.72; vsize = [1 1 1 TR]; vtype = 'f'; filename = 'simple';
data = struct;
for s = 1:S
    for r = 1
        EC = reshape(Y((s-1)*t+1:s*t,:),10,10,V/100,t);
        save_avw(EC,sprintf('Results/Temp/sub%02d_run%02d_%s',s,r,filename),vtype,vsize);
        data.(sprintf('S%02d',s)).(sprintf('R%02d',r)) = sprintf('%s/Results/Temp/sub%02d_run%02d_%s.nii.gz',pwd,s,r,filename);
    end
end
saveJSONfile(data,'Results/simple.json');

% Run PROFUMO
system('nice -n 20 ~samh/bin/PROFUMO Results/simple.json 2 Results/PROFUMO_simple --useHRF 0.72 --hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf -d 0.05 > Results/Output_PROFUMO_simple.txt')

% Load PROFUMO results
Pdir = 'Results/PROFUMO_simple.pfm/FinalModel/';
PFMgmap = PFM_loadGroupSpatialMaps_new(Pdir);
PFMmaps = cell(S,1); PFMts = cell(S,1); PFMnet = cell(S,1);
for s = 1:S
    means = h5read(sprintf('%s/Subjects/S%02d/SpatialMaps.post/Signal/Means.hdf5',Pdir,s),'/dataset');
    probs = h5read(sprintf('%s/Subjects/S%02d/SpatialMaps.post/Signal/MembershipProbabilities.hdf5',Pdir,s),'/dataset');
    PFMmaps{s} = means .* probs; clear means probs
    PFMnet{s} = h5read(sprintf('%s/Subjects/S%02d/TemporalPrecisionMatrix.post/Mean.hdf5',Pdir,s),'/dataset');
    PFMts{s} = demean(h5read(sprintf('%s/Subjects/S%02d/Runs/R01/TimeCourses.post/CleanTimeCourses/Means.hdf5',Pdir,s),'/dataset'));
end 
