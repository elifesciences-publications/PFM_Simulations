function saveNIfTIs(D, params, directory)

addpath('IO/saveJSONfile/')

% Prepare nifti
vsize = [1 1 1 params.TR];
vtype = 'f';

% Loop over subjects to save nifti & create data structure of json save
data = struct;
%SubjectList = struct;
%RunList = struct; I = 1;
for s = 1:params.S
    subj = sprintf('S%02d',s);
    
    %SubjectList.(subj) = [];
    
    for r = 1:params.R(1)
        run = sprintf('R%02d',r);
        
        data3d = reshape(D{s}{r}, 10, 10, params.V / 100, params.T);
        fileName = fullfile(pwd(), directory, [subj '_' run '.nii.gz']);
        save_avw(data3d, fileName, vtype, vsize);
        data.(subj).(run) = fileName;
        
        %RunList.(sprintf('i%02d',I)).Subject = subj;
        %RunList.(sprintf('i%02d',I)).Run = run; I = I+1;
    end
end

% Save json
saveJSONfile(data, fullfile(directory, 'PROFUMO_SpecFile.json'))
%saveJSONfile(SubjectList, fullfile(directory, 'PROFUMO_SubjectList.json'))
%saveJSONfile(RunList, fullfile(directory, 'PROFUMO_RunList.json'))

return
