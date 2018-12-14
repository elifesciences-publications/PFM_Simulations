% Generates a test fMRI data set and compares the performance of a variety of node discovery algorithms
% Compares results of PROFUMO and ICA
% Main aim is to test effects of spatial overlap vs spatial misalignment

close all; clear all; clc

% Inputs
fileName = 'Simulations';
overlap = true
misalignment = true

% Fix filename
if ~overlap
    fileName = strcat(fileName, '_NoOverlap');
end
if ~misalignment
    fileName = strcat(fileName, '_NoMisalignment');
end
fileName

% Plotting
plotFigures = false;

%% Set paths

restoredefaultpath

addpath(genpath('~samh/Documents/Code/MATLAB/'))
addpath(strcat(getenv('FSLDIR'), '/etc/matlab/'))

% Internal paths
addpath('DataGeneration/');
addpath('Methods/');
addpath('Scoring/');
addpath('Visualisation/');
addpath('Overlap_functions/');

if plotFigures
    prettyFigures();
end

%% Set size of problem

%Number of times to repeat simulation / test cycle
params.nRepeats = 10;

%Details of scans
params.S = 30;       %Subjects
params.R = 2*ones(params.S,1);   %Runs

params.T = 600;     %No. of time points per fMRI scan
params.TR = 0.72;
params.dt = 0.1;     %Neural sampling rate
%Amount of neural points to simulate - more than scan length so don't have
%to zero pad HRF at start of scan
params.Tn = ceil(1.25 * params.T * params.TR / params.dt);

%% Set size of atlas / mode matrix

% Atlas
atlasParams = params;
atlasParams.V = 10000;    %Voxels
atlasParams.N = 100;      %Number of nodes in the atlas

% Modes
modeParams = params;
modeParams.V = atlasParams.N;
modeParams.N = 15;        %Number of modes

% And store appropriate values for combined maps
params.N = modeParams.N;
params.V = atlasParams.V;

%Number of modes to infer
params.iN = 15;

%% Set the details of the tests we want to do

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Atlas

atlasOptions.P.form = 'Probabilistic';

%Choose form for Pg
atlasOptions.Pg.form = 'BlockAtlas';
% widthPrecision controls the variability in the size of the parcels
% Smaller numbers give more variability in parcel sizes
atlasOptions.Pg.widthPrecision = 25;
% Post-hoc smoothing of maps (width, in voxels, of the filter)
atlasOptions.Pg.smootherWidth = 0.1 * (atlasParams.V / atlasParams.N);

%Choose form for Ps
atlasOptions.Ps.form = 'WeightedGamma';
% Probability subject voxel is not drawn from group distribution
if overlap
    atlasOptions.Ps.p = 0.0005;
else
    atlasOptions.Ps.p = 0.0;
end
% Minimum weight - useful to make sure all weights are different from noise
atlasOptions.Ps.minWeight = 0.1;
% Weights are gamma(a,b) distributed (mean = a/b)
% Increasing a,b makes them more stable
atlasOptions.Ps.weightRange.a = 0.9 * 20.0;
atlasOptions.Ps.weightRange.b = 20.0;
% Little bit of Gaussian noise for old times sake
atlasOptions.Ps.epsilon = 0.025;

% Choose registration
if misalignment
    atlasOptions.P.registration.form = 'RandomSmooth';
    atlasOptions.P.registration.maxError = 1.5 * (atlasParams.V / atlasParams.N);
    % This parameter controls the size of misalignments
    % It represents the furthest one voxel can be moved by misregistration
    % Useful to express this in terms of `c * (atlas.V / atlas.N)`, i.e. the
    % average parcel size. The parameter `c` then represents the max
    % misalignment in terms of number of parcels rather than voxels.
else
    atlasOptions.P.registration.form = 'Null';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modes

modeOptions.P.form = 'Probabilistic';

modeOptions.Pg.form = 'BiasedBoxcar';
% How many spatially contiguous blocks per mode? Follows `Poisson(nBlocks) + 1`
modeOptions.Pg.nBlocks = 0.5;
% How big are the modes? On average, they cover `p * V` voxels
% If we have N modes, then we expect `p * N` modes in every voxel
% This is therefore a crude proxy for overlap
%%% HighOverlap: 1.4; LowOverlap 1.2; %%%
if overlap
    modeOptions.Pg.p = 1.3 / params.N;
else
    modeOptions.Pg.p = 0.95 / params.N;
end
modeOptions.Pg.pVar = 0.01 ^ 2; % i.e. p will vary over approximately +/- 2.0 * sqrt(pVar)
% Proportion of (secondary) blocks that are positive
modeOptions.Pg.pPosBlock = 0.5;

%Choose form for Ps
modeOptions.Ps.form = 'WeightedGamma';
% Probability subject voxel is not drawn from group distribution
% `p = c / (V * N)` means that, on average `c` parcels are active
% in a given subject that were not in the group maps
if overlap
    modeOptions.Ps.p = 2.0 / (modeParams.V * modeParams.N);
else
    modeOptions.Ps.p = 0.0;
end
% Minimum weight - useful to make sure all weights are different from noise
modeOptions.Ps.minWeight = 0.0;
% Weights are gamma(a,b) distributed (mean = a/b)
% Increasing a,b makes them more stable
modeOptions.Ps.weightRange.a = 5.0;
modeOptions.Ps.weightRange.b = 5.0;
% Little bit of Gaussian noise for old times sake
modeOptions.Ps.epsilon = 0.01;

%Choose registration errors
modeOptions.P.registration.form = 'Null';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time courses
options.An.form = 'Freq';
options.An.offRate = 0.05 * 1/params.N;
switch options.An.form
    case 'Freq'
        % Reducing these parameters will increase the
        % strength of the correlations at the group,
        % subject and run level respectively
        options.An.Cg_dof = 250;
        options.An.Cs_dof = 500;
        options.An.Cr_dof = 500;
        options.An.p = 0.1;
        options.An.fc = 0.1; %in Hz
        options.An.fAmp = 2;
        options.An.epsilon = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOLD signal
options.BS.form = 'SaturatingFlobsHRF';
switch options.BS.form
    case 'Linear'
        
    case 'FlobsHRF'
        options.BS.HRFcoeffs.mu = [1 0 -0.2];
        options.BS.HRFcoeffs.sigma = [0.1 0.1 0.1];
        options.BS.nongaussian = 0;
        
    case 'SaturatingFlobsHRF'
        options.BS.HRFcoeffs.mu = [1 0 -0.2];
        options.BS.HRFcoeffs.sigma = [0.1 0.1 0.1];
        options.BS.tanhPercentile = 99;
        options.BS.tanhMax = 0.9;
        options.BS.nongaussian = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noise

%Signal to noise ratio (expressed in terms of power)
options.D.SNR = 0.1;

%Type
options.D.noise.form = 'StructuredTdist';
switch options.D.noise.form
    case 'SpatiotemporallyWhiteGaussian'
        
    case 'SpatiotemporallyWhiteTdist'
        options.D.noise.a = 8;
        options.D.noise.b = 8;
        
    case 'StructuredTdist'
        options.D.noise.a = 3.0;
        options.D.noise.b = 3.0;
        % Standard deviations of the respective subspaces
        options.D.noise.structuredStd = 1.0;
        options.D.noise.unstructuredStd = 4.0;
        % Rank of structured noise
        options.D.noise.N = 5;
end


%% Run the tests

rng('shuffle')

for n = 1:params.nRepeats
    
    fileNameN = sprintf('%s_%02d', fileName, n)
    
    %% Generate data
    
    if plotFigures && (n==1)
        plotNow = true;
    else
        plotNow = false;
    end
    
    atlasP = generateMaps(atlasParams, atlasOptions, plotNow);
    
    [modeP, plotMaps] = generateMaps(modeParams, modeOptions, plotNow);
    
    %Combine to make full maps
    P = cell(params.S,1);
    for s = 1:params.S
        P{s} = atlasP{s} * modeP{s};
    end
    %Plot if requested
    if plotNow
        plotMaps(P, params, options);
    end
    
%     %Combine to make full maps
%     addpath('~steve/matlab/icasso122','~steve/matlab/FastICA_25');
%     addpath /vols/Scratch/janineb/HCP/DMN1200/Functions
%     P = cell(params.S,1);
%     for s = 1:params.S
%         P{s} = atlasP{s} * modeP{s};
%         [icaS,icaA,icaW] = fastica(P{s}','approach','symm','g','tanh','epsilon',1e-11,'maxNumIterations',3000,'lastEig',15);
%         [C12,munkres_assign] = spatialcorr(icaS',P{1});
%         [i,j] = find(munkres_assign==1); C12 = sign(C12(munkres_assign==1));
%         icaS = icaS(i,:)'; icaS = icaS.*repmat(C12',size(icaS,1),1);
%         P{s} = icaS;
%     end
    
    An = generateNeuralTimecourses(params, options, plotNow);
    
    PA = generateBoldSignal(P, An, params, options, plotNow);
    
    D = generateData(PA, params, options, plotNow);
    
    %Finally, add a global rescaling such that all scans are overall
    %unit variance
    vD = 0;
    for s = 1:params.S
        for r = 1:params.R(s)
            vD = vD + var(D{s}{r}(:));
        end
    end
    vD = vD / sum(params.R);
    for s = 1:params.S
        for r = 1:params.R(s)
            D{s}{r} = D{s}{r} / sqrt(vD);
            PA{s}{r} = PA{s}{r} / sqrt(vD);
        end
    end
    
    %% Look at ground truth accuracy in the PA subspace
    % This looks at how well the linear mixing model can do, in the best
    % case scenario
    
    %Extract mean group map
    Pg = 0;
    for s = 1:params.S
        Pg = Pg + P{s};
    end
    Pg = Pg / params.S;
    
    %Use the group and the subject specific maps to regress the timecourses
    %out of the BOLD signal (noise free)
    A = cell(params.S,1); Ag = cell(params.S,1);
    for s = 1:params.S
        A{s} = cell(params.R(s),1); Ag{s} = cell(params.R(s),1);
        for r = 1:params.R(s)
            Ag{s}{r} = Pg \ PA{s}{r};
            A{s}{r} = P{s} \ PA{s}{r};
        end
    end
    
    %Save scores
    %Ground truth linear model
    scores.GT.PA(n) = calculateBoldRecovery(PA, makePA(P,A,params), params);
    
    %Ground truth mean map model
    scores.GTg.PA(n) = calculateBoldRecovery(PA, ...
        makePA(repmat({Pg},params.S,1),Ag,params), params);
    
    %Mean map accuracy (if we are trying to find the true number of modes)
    if params.N == params.iN
        [scores.GTg.P(:,n), scores.GTg.A(:,n), scores.GTg.cP(:,n), ...
            scores.GTg.cA(:,n)] = calculateDecompositionAccuracy( ...
            P, repmat({Pg}, params.S, 1), A, Ag, params);
    end
    
    %% Save simulated data
    save(sprintf('Results/PFMsims_atlas_%s',fileNameN), ...
        'P', 'Pg', 'A', 'D', ...
        'atlasParams', 'modeParams', 'params', 'scores', ...
        '-v7.3')
    
    %% Save NIFTIs and run latest versions of MELODIC and PROFUMO
    
    Save_niftis(D,params,fileNameN,atlasParams);
    system(sprintf('sh Overlap_functions/ICA_PROFUMO.sh %1.2f %s %d %s',params.TR,fileNameN,params.iN,pwd))
    
end

%% Plot results

if plotFigures
    plotScores(scores, params);
    input('Press return to continue')
    close all
end

% When finished - run next script to produce figure for paper
Paper_figure_new({fileName});
