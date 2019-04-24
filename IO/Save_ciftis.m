function Save_ciftis(D,params,filename)

% Set paths
path(path,'/home/fs0/janineb/scratch/matlab/cifti-matlab/')
path(path,'/home/fs0/janineb/scratch/matlab/saveJSONfile/')

% Load Data
load(sprintf('Results/PFMsims_atlas_%s.mat',filename))
S = size(D{1}{1});
TR = params.TR;

% Prepare cifti
EC = ft_read_cifti('../PFM/maps_full_originalOrder.dtseries.nii');
EC.hdr.dim(7) = S(1); EC.hdr.dim(6) = S(2);
EC.time = 1:TR:params.T*TR+0.5;
EC.pos = EC.pos(1:S(1),:);
EC = rmfield(EC,{'brainstructure','brainstructurelabel','dim','transform'});

% Loop over subjects to save cifti & create data structure of json save
data = struct;
SubjectList = struct;
RunList = struct; I = 1;
for s = 1:params.S
    SubjectList.(sprintf('S%02d',s)) = [];
    for r = 1:params.R(1)
        EC.dtseries = D{s}{r};
        ft_write_cifti(sprintf('Results/Temp/sub%02d_run%02d_%s',s,r,filename),EC,'parameter','dtseries');
        data.(sprintf('S%02d',s)).(sprintf('R%02d',r)) = sprintf('/vols/Scratch/janineb/HCP/DMN1200/PFM_simulations/Results/Temp/sub%02d_run%02d_%s.dtseries.nii',s,r,filename);
        RunList.(sprintf('i%02d',I)).Subject = sprintf('S%02d',s);
        RunList.(sprintf('i%02d',I)).Run = sprintf('R%02d',r); I = I+1;
    end
end

% Save json
saveJSONfile(data,sprintf('Results/PFMsims_atlas_%s.json',filename))
% saveJSONfile(SubjectList,sprintf('Results/PFMsims_atlas_%s_SubjectList.json',Option))
% saveJSONfile(RunList,sprintf('Results/RunList_S%02d_R%02d.json',params.S,params.R(1)))

