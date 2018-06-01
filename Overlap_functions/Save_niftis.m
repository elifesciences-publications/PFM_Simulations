function Save_niftis(D,params,filename,atlasParams)

addpath('Overlap_functions/Dependencies/saveJSONfile/')

% Load Data
S = size(D{1}{1});
TR = params.TR;

% Prepare nifti
vsize = [1 1 1 TR];
vtype = 'f';

% Loop over subjects to save cifti & create data structure of json save
data = struct;
SubjectList = struct;
RunList = struct; I = 1;
for s = 1:params.S
    SubjectList.(sprintf('S%02d',s)) = [];
    for r = 1:params.R(1)
        EC = reshape(D{s}{r},10,10,atlasParams.V/100,params.T);
        save_avw(EC,sprintf('Results/Temp/sub%02d_run%02d_%s',s,r,filename),vtype,vsize);
        data.(sprintf('S%02d',s)).(sprintf('R%02d',r)) = sprintf('%s/Results/Temp/sub%02d_run%02d_%s.nii.gz',pwd,s,r,filename);
        RunList.(sprintf('i%02d',I)).Subject = sprintf('S%02d',s);
        RunList.(sprintf('i%02d',I)).Run = sprintf('R%02d',r); I = I+1;
    end
end

% Save json
saveJSONfile(data,sprintf('Results/PFMsims_atlas_%s.json',filename))
% saveJSONfile(SubjectList,sprintf('Results/PFMsims_atlas_%s_SubjectList.json',Option))
% saveJSONfile(RunList,sprintf('Results/RunList_S%02d_R%02d.json',params.S,params.R(1)))

