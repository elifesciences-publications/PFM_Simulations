% Generates a test fMRI data set and compares the performance of a variety of node discovery algorithms
% Compares results of PROFUMO and ICA
% Main aim is to test effects of spatial overlap vs spatial misalignment

close all; clear all; clc

% Inputs
fileName = 'test';
No_overlap = true; 
plot_SamFigures = false;
plot_JanineFigures = true;

%% Set paths
restoredefaultpath
addpath(genpath('~samh/Documents/Code/Simulated_fMRI_Tests/'))
addpath(genpath('~samh/Documents/Code/MATLAB/'))
addpath(genpath('~samh/Documents/Code/Algorithms/'))
addpath('~samh/Documents/Code/VBGP_Tests/')
addpath('DataGeneration/'); 
addpath('Methods/');
addpath('Scoring/'); 
addpath('Visualisation/')
addpath('Overlap_functions/');
addpath('/opt/fmrib/fsl/etc/matlab/')

rng('shuffle')
if plot_SamFigures
    prettyFigures()
end

%% Set size of problem

%Number of times to repeat simulation / test cycle
params.nRepeats = 1;

%Details of scans
params.S = 30;       %Subjects
params.R = 2*ones(params.S,1);   %Repeats

params.T = 600;     %No. of time points per fMRI scan
params.TR = 0.72;
params.dt = 0.1;     %Neural sampling rate
%Amount of neural points to simulate - more than scan length so don't have
%to zero pad HRF at start of scan
params.Tn = ceil(1.25 * params.T * params.TR / params.dt);

%% Set size of atlas / mode matrix

atlasParams = params;
modeParams = params;

atlasParams.V = 2500;     %Voxels
atlasParams.N = 100;      %Number of nodes in the atlas

modeParams.V = atlasParams.N;
modeParams.N = 15;        %Number of modes

params.N = modeParams.N;
params.V = atlasParams.V;

%Number of modes to infer
params.iN = 15;
%params.iN = 25
%params.iN = 40

%% Set the details of the tests we want to do

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Atlas

atlasOptions.P.form = 'Additive';

%Choose form for Pg
atlasOptions.Pg.form = 'BlockAtlas';
% widthPrecision controls the variability in the size of the parcels
% Smaller numbers give more variability in parcel sizes
atlasOptions.Pg.widthPrecision = 25;

%Choose form for Ps
%atlasOptions.Ps.form = 'Null';
atlasOptions.Ps.form = 'Gaussian';
atlasOptions.P.PsPg = 0.10; %Ratio of std(P{s}(:)-Pg(:)) / std(Pg(:))

% Choose registration
atlasOptions.P.registration.form = 'RandomSmooth';
atlasOptions.P.registration.maxError = 1.5 * (atlasParams.V / atlasParams.N);
% This parameter controls the size of misalignments
% It represents the furthest one voxel can be moved by misregistration
% Useful to express this in terms of `c * (atlas.V / atlas.N)`, i.e. the
% average parcel size. The parameter `c` then represents the max
% misalignment in terms of number of parcels rather than voxels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modes

modeOptions.P.form = 'Additive';

%Choose form for Pg
if No_overlap
    modeOptions.Pg.form = 'BlockAtlas';
    % widthPrecision controls the variability in the size of the parcels
    % Smaller numbers give more variability in parcel sizes
    modeOptions.Pg.widthPrecision = 10; %25;
else
    modeOptions.Pg.form = 'BiasedBoxcar';
    % How many spatially contiguous blocks per mode? Follows Poisson(nBlocks)
    modeOptions.P.nBlocks = 2; %4;
    % How big are the modes? On average, they cover `p * V` voxels
    % If we have N modes, then we expect `p * N` modes in every voxel
    % This is therefore a crude proxy for overlap
    modeOptions.P.p = 1.0 / params.N;
    modeOptions.P.pVar = 0.00075;
    % Increase this parameter to make blocks less likely to overlap
    modeOptions.P.biasStrength = 0.75;
    % Minimum weight - useful to make sure all weights are different from noise
    modeOptions.P.minWeight = 0.5;
    % Weights are gamma(a,b) distributed (mean = a/b)
    % Increasing a,b makes them more stable
    modeOptions.P.weightRange.a = 5;
    modeOptions.P.weightRange.b = 5;
    % Proportion of blocks that are positive
    modeOptions.P.pPosBlock = 0.8;
    % Post-hoc smoothing of maps
    modeOptions.P.smootherWidth = floor( 0.01*modeParams.V );
end

%Choose form for Ps
modeOptions.Ps.form = 'SS';
modeOptions.P.PsPg = 0.05; %Ratio of std(P{s}(:)-Pg(:)) / std(Pg(:))
%Set options based on that choice
switch modeOptions.Ps.form
    
    case 'SS'
        % Strength of correlations in the noise
        modeOptions.Ps.rot = 0.25;
        % Noise sparsity
        modeOptions.Ps.p = 0.05;
        % Standard deviations of spike and slab
        modeOptions.Ps.sigma = 1;
        modeOptions.Ps.epsilon = 0.1;
        
    case 'Gaussian'
        
    case 'Null'
end

%Choose registration errors
modeOptions.P.registration.form = 'Null';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time courses
options.An.form = 'Freq';
options.An.offRate = 0.25 * 1/params.N;
switch options.An.form
    case 'Freq';
        % Increasing these parameters will increase the
        % strength of the correlations at the group,
        % subject and run level respectively
        options.An.rot = 0.3;
        options.An.rotS = 0.5;
        options.An.rotR = 0.1;
        options.An.p = 0.2;
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
        options.BS.HRFcoeffs.mu = [1 0 0];
        options.BS.HRFcoeffs.sigma = [0.1 0.1 0.1];
        
    case 'SaturatingFlobsHRF'
        options.BS.HRFcoeffs.mu = [1 0 0];
        options.BS.HRFcoeffs.sigma = [0.1 0.1 0.1];
        options.BS.tanhPercentile = 99;
        options.BS.tanhMax = 0.9;
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
        options.D.noise.a = 8;
        options.D.noise.b = 8;
        % Standard deviations of the respective subspaces
        options.D.noise.structuredStd = 0.25;
        options.D.noise.unstructuredStd = 0.75;
        % Rank of structured noise
        options.D.noise.N = 20;
end


%% Run the tests

for n = 1:params.nRepeats
    
    %% Generate data
    
    if plot_SamFigures && (n==1)
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
    
    if plot_SamFigures
        cP = P{1}' * P{1};
        cP = cP ./ sqrt(diag(cP)*diag(cP)');
        figure; imagesc(cP, [-1 1]); colorbar; colormap(bluewhitered)
        xlabel('Mode'); ylabel('Mode');
        axis square
        set(gcf, 'Position', [200 200 500 400])
        
        cA = [An{1}{1} An{1}{2} An{1}{3} An{1}{4}] * [An{1}{1} An{1}{2} An{1}{3} An{1}{4}]';
        cA = cA ./ sqrt(diag(cA)*diag(cA)');
        figure; imagesc(cA, [-1 1]); colorbar; colormap(bluewhitered)
        xlabel('Mode'); ylabel('Mode');
        axis square
        set(gcf, 'Position', [200 200 500 400])
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
    
    %% Extract PFMs
    
    %For the first run, may want to plot convergence
    if plotNow
        %If we want to run with plots of convergence, pass in best guesses
        %for mean map and node time courses
        [pfmP1, pfmA1, pfmPg1] = runVBGP(D, params, params.iN, Pg, A);
    else
        %Just run VBGP normally
        [pfmP1, pfmA1, pfmPg1] = runVBGP(D, params, params.iN);
    end
    
    
    %Save scores
    scores.PFMs.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(pfmP1,pfmA1,params), params);
    [scores.PFMs.P(:,2*(n-1)+1), scores.PFMs.A(:,2*(n-1)+1), ...
        scores.PFMs.cP(:,2*(n-1)+1), scores.PFMs.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, pfmP1, A, pfmA1, params);
    
    
    %% Run PCA (using SVD)
    
    [pcaP, pcaA, svdU, svdS, svdV] = runPCA(D, params, params.iN);
    
    %fprintf('Hello, Sam!! \n');
    
    %Save scores
    scores.PCA.PA(n) = calculateBoldRecovery(PA, ...
        makePA(repmat({pcaP},params.S,1),pcaA,params), params);
    
    [scores.PCA.P(:,n), scores.PCA.A(:,n), scores.PCA.cP(:,n), ...
        scores.PCA.cA(:,n)] = calculateDecompositionAccuracy( ...
        P, repmat({pcaP}, params.S, 1), A, pcaA, params);
    
    %% Run sICA on low dim PCA, with dual reg
    
    [sicaP1, sicaA1] = runSICA(svdU, svdS, svdV, params, params.iN);
    sicaP1 = repmat({sicaP1},params.S,1);
    
    %Save scores
    scores.sICA.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(sicaP1,sicaA1,params), params);
    [scores.sICA.P(:,2*(n-1)+1), scores.sICA.A(:,2*(n-1)+1), ...
        scores.sICA.cP(:,2*(n-1)+1), scores.sICA.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, sicaP1, A, sicaA1, params);
    
    %Run dual regression
    [sicaP1_DR, sicaA1_DR] = runDR(D, sicaA1, params);
    
    %%% ADDED BY JANINE: TRADITIONAL DR
    addpath /usr/local/fmrib/fmt/
    sica_P1_DRtest = cell(params.S,1);
    sica_A1_DRtest = cell(params.S,1);
    for s = 1:params.S
        sica_A1_DRtest{s} = cell(params.R(s),1);
        Ds = zeros(params.V, params.T*params.R(s));
        M = zeros(atlasParams.V,params.iN,params.R(1));
        for r = 1:params.R(s)
            Ds(:, (r-1)*params.T+(1:params.T)) = D{s}{r};
            sica_A1_DRtest{s}{r} = demean((pinv(demean(double(sicaP1{1})))*demean(double(D{s}{r})))');
            M(:,:,r) = demean((pinv(demean(double(sica_A1_DRtest{s}{r})))*demean(double(D{s}{r}))')');
        end
        sica_P1_DRtest{s} = mean(M,3);
    end
    
    %Save scores
    scores.sICA_DR.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(sicaP1_DR,sicaA1_DR,params), params);
    [scores.sICA_DR.P(:,2*(n-1)+1), scores.sICA_DR.A(:,2*(n-1)+1), ...
        scores.sICA_DR.cP(:,2*(n-1)+1), scores.sICA_DR.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, sicaP1_DR, A, sicaA1_DR, params);
      
    %% Run tICA on low dim PCA, with dual reg
    
    [ticaP1, ticaA1] = runTICA(svdU, svdS, svdV, params, params.iN);
    ticaP1 = repmat({ticaP1},params.S,1);
    
    %Save scores
    scores.tICA.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(ticaP1,ticaA1,params), params);
    [scores.tICA.P(:,2*(n-1)+1), scores.tICA.A(:,2*(n-1)+1), ...
        scores.tICA.cP(:,2*(n-1)+1), scores.tICA.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, ticaP1, A, ticaA1, params);
    
    %Run dual regression
    [ticaP1_DR, ticaA1_DR] = runDR(D, ticaA1, params);
    
    %Save scores
    scores.tICA_DR.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(ticaP1_DR,ticaA1_DR,params), params);
    [scores.tICA_DR.P(:,2*(n-1)+1), scores.tICA_DR.A(:,2*(n-1)+1), ...
        scores.tICA_DR.cP(:,2*(n-1)+1), scores.tICA_DR.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, ticaP1_DR, A, ticaA1_DR, params);
       
    %% Run sICA & tICA
    % sICA on high dim PCA, with dual reg, followed by tICA
    
    [sticaP1, sticaA1] = runSTICA(D, svdU, svdS, svdV, params, atlasParams.N, params.iN);
    
    %Save scores
    scores.stICA.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(sticaP1,sticaA1,params), params);
    [scores.stICA.P(:,2*(n-1)+1), scores.stICA.A(:,2*(n-1)+1), ...
        scores.stICA.cP(:,2*(n-1)+1), scores.stICA.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, sticaP1, A, sticaA1, params);
    
    %Run dual regression
    [sticaP1_DR, sticaA1_DR] = runDR(D, sticaA1, params);
    
    %Save scores
    scores.stICA_DR.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(sticaP1_DR,sticaA1_DR,params), params);
    [scores.stICA_DR.P(:,2*(n-1)+1), scores.stICA_DR.A(:,2*(n-1)+1), ...
        scores.stICA_DR.cP(:,2*(n-1)+1), scores.stICA_DR.cA(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, sticaP1_DR, A, sticaA1_DR, params);
    
    %% Run postprocessing 
    Save_niftis(D,params,fileName,atlasParams);
    system(sprintf('sh Overlap_functions/ICA_PROFUMO.sh %1.2f %s %d %s',params.TR,fileName,params.iN,pwd))
    %unix(sprintf('nice -n 20 ~samh/bin/PROFUMO Results/PFMsims_atlas_%s.json %d Results/PROFUMO_PFMsim_atlas_%s --useHRF %1.2f --hrfFile ~samh/PROFUMO/Scripts/DefaultHRF.phrf -d 0.05 > Results/Output_PROFUMO_PFMsims_atlas_%s.txt',fileName,params.iN,fileName,params.TR,fileName))
    pause(300)
    [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new] = Melodic_DR(fileName,D,atlasParams,params);
    [pfmPg1_new,pfmP1_new,pfmA1_new] = loadNewPROFUMO(fileName,params);
    if plot_JanineFigures
        RunOverlapTests_plots
    end
    
    %% Save results
    save(sprintf('Results/PFMsims_atlas_%s',fileName),'-v7.3')
    
end

%% Plot results

if plot_SamFigures
    plotScores(scores, params);
    input('Press return to continue')
    close all
end

