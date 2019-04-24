function [pfmPg1_new,pfmP1_new,pfmA1_new,pfmNet] = loadNewPROFUMO(filename,params) 

% Set paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/saveJSONfile/')
path(path,'/vols/Scratch/janineb/HCP/DMN1200/Functions')

% Load Data
Pdir = sprintf('Results/PROFUMO_PFMsims_atlas_%s.pfm/FinalModel/',filename);

pfmPg1_new = PFM_loadGroupSpatialMaps_new(Pdir);

% Replace with New PROFUMO results
pfmP1_new = cell(params.S,1); pfmA1_new = cell(params.S,1);
pfmNet = cell(params.S,1);
for s = 1:params.S
    means = h5read(sprintf('%s/Subjects/S%02d/SpatialMaps.post/Signal/Means.hdf5',Pdir,s),'/dataset');
    probs = h5read(sprintf('%s/Subjects/S%02d/SpatialMaps.post/Signal/MembershipProbabilities.hdf5',Pdir,s),'/dataset');
    pfmP1_new{s} = means .* probs; clear means probs
    pfmNet{s} = h5read(sprintf('%s/Subjects/S%02d/TemporalPrecisionMatrix.post/Mean.hdf5',Pdir,s),'/dataset');
    for r = 1:params.R(1)
        pfmA1_new{s}{r} = demean(h5read(sprintf('%s/Subjects/S%02d/Runs/R%02d/TimeCourses.post/CleanTimeCourses/Means.hdf5',Pdir,s,r),'/dataset'));
    end
end
    


