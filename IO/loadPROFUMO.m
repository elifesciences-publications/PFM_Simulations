function [ Pg, P, A ] = loadPROFUMO( baseFileName, params )

pfmDir = fullfile([baseFileName '_PROFUMO.pfm'], 'FinalModel');

% Load group maps
mPg = h5read(fullfile(pfmDir, 'GroupSpatialModel', 'SignalMeans.post', ...
    'Means.hdf5'), '/dataset');
pPg = h5read(fullfile(pfmDir, 'GroupSpatialModel', 'Memberships.post', ...
    'Class_2', 'Probabilities.hdf5'), '/dataset');
Pg = mPg .* pPg; clear mPg pPg;

% And subject maps / timecourses
P = cell(params.S, 1);
A = cell(params.S, 1);
%netmats = cell(params.S,1);
for s = 1:params.S
    subj = sprintf('S%02d',s);
    subjDir = fullfile(pfmDir, 'Subjects', subj);
    
    % Subject maps
    mPs = h5read(fullfile(subjDir, 'SpatialMaps.post', 'Signal', ...
        'Means.hdf5'), '/dataset');
    pPs = h5read(fullfile(subjDir, 'SpatialMaps.post', 'Signal', ...
        'MembershipProbabilities.hdf5'), '/dataset');
    P{s} = mPs .* pPs; clear mPs pPs;
    
    % Time courses
    A{s} = cell(params.R(s), 1);
    for r = 1:params.R(s)
        run = sprintf('R%02d',r);
        runDir = fullfile(subjDir, 'Runs', run);
        
        A{s}{r} = h5read(fullfile(runDir, 'TimeCourses.post', ...
            'CleanTimeCourses', 'Means.hdf5'), '/dataset');
    end
    
    % Temporal netmats
    %netmats{s} = h5read(fullfile(subjDir, 'TemporalPrecisionMatrix.post', ...
    %    'Mean.hdf5'), '/dataset');
end
