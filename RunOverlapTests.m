% Generates a test fMRI data set and compares the performance of a variety of node discovery algorithms
% Compares results of PROFUMO and ICA
% Main aim is to test effects of spatial overlap vs spatial misalignment

close all; clear all; clc

% Inputs
fileName = 'LowO_LowM_MedT_NoHighs';
No_overlap = false;
plot_SamFigures = false;
plot_JanineFigures = false;

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

atlasParams = params;
modeParams = params;

atlasParams.V = 10000;     %Voxels
atlasParams.N = 100;      %Number of nodes in the atlas

modeParams.V = atlasParams.N;
modeParams.N = 25;        %Number of modes

params.N = modeParams.N;
params.V = atlasParams.V;

%Number of modes to infer
params.iN = 15;
%params.iN = 25
%params.iN = 40

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
atlasOptions.Pg.smootherWidth = 25.0;

%Choose form for Ps
atlasOptions.Ps.form = 'WeightedGamma';
% Probability subject voxel is not drawn from group distribution
atlasOptions.Ps.p = 0.001;
% Minimum weight - useful to make sure all weights are different from noise
atlasOptions.Ps.minWeight = 0.1;
% Weights are gamma(a,b) distributed (mean = a/b)
% Increasing a,b makes them more stable
atlasOptions.Ps.weightRange.a = 0.9 * 20.0;
atlasOptions.Ps.weightRange.b = 20.0;
% Little bit of Gaussian noise for old times sake
atlasOptions.Ps.epsilon = 0.01;

% Choose registration
atlasOptions.P.registration.form = 'RandomSmooth';
atlasOptions.P.registration.maxError = 0.75 * (atlasParams.V / atlasParams.N);
% This parameter controls the size of misalignments
% It represents the furthest one voxel can be moved by misregistration
% Useful to express this in terms of `c * (atlas.V / atlas.N)`, i.e. the
% average parcel size. The parameter `c` then represents the max
% misalignment in terms of number of parcels rather than voxels.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modes

modeOptions.P.form = 'Probabilistic';

modeOptions.Pg.form = 'BiasedBoxcar';
% How many spatially contiguous blocks per mode? Follows `Poisson(nBlocks) + 1`
modeOptions.Pg.nBlocks = 1.25;
% How big are the modes? On average, they cover `p * V` voxels
% If we have N modes, then we expect `p * N` modes in every voxel
% This is therefore a crude proxy for overlap
%%% HighOverlap: 1.4; LowOverlap 1.2; %%%
modeOptions.Pg.p = 1.2 / params.N;
modeOptions.Pg.pVar = 0.01 ^ 2; % i.e. p will vary over approximately +/- 2.0 * sqrt(pVar)
% Increase this parameter to make blocks less likely to overlap
% Between 0 and 1
%%% HighOverlap: 0.5; LowOverlap 0.9; %%%
modeOptions.Pg.biasStrength = 0.75;
% Proportion of (secondary) blocks that are positive
modeOptions.Pg.pPosBlock = 0.7;

%Choose form for Ps
modeOptions.Ps.form = 'WeightedGamma';
% Probability subject voxel is not drawn from group distribution
% `p = c / (V * N)` means that, on average `c` parcels are active
% in a given subject that were not in the group maps
modeOptions.Ps.p = 5.0 / (modeParams.V * modeParams.N);
% Minimum weight - useful to make sure all weights are different from noise
modeOptions.Ps.minWeight = 0.0;
% Weights are gamma(a,b) distributed (mean = a/b)
% Increasing a,b makes them more stable
modeOptions.Ps.weightRange.a = 2.0;
modeOptions.Ps.weightRange.b = 2.0;
% Little bit of Gaussian noise for old times sake
modeOptions.Ps.epsilon = 0.01;

%Choose registration errors
modeOptions.P.registration.form = 'Null';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time courses
options.An.form = 'Freq';
options.An.offRate = 0.05 * 1/params.N;
switch options.An.form
    case 'Freq';
        % Increasing these parameters will increase the
        % strength of the correlations at the group,
        % subject and run level respectively
        options.An.rot = 0.3; % 0.3
        options.An.rotS = 0.5; % 0.5
        options.An.rotR = 0.1; % 0.1
        options.An.p = 0.2;
        options.An.fc = 0.1; %in Hz
        options.An.fAmp = 2;
        options.An.epsilon = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BOLD signal
options.BS.form = 'FlobsHRF';
switch options.BS.form
    case 'Linear'
        
    case 'FlobsHRF'
        options.BS.HRFcoeffs.mu = [1 0 0];
        options.BS.HRFcoeffs.sigma = [0.1 0.1 0.1];
        options.BS.nongaussian = 1;
        
    case 'SaturatingFlobsHRF'
        options.BS.HRFcoeffs.mu = [1 0 0];
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
n=1;
while n <= params.nRepeats
    
    fileNameN = sprintf('%s_%02d',fileName,n);
    fprintf('Trying\n');
    
    %% Generate data
    params.N = 25;
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
    
    % Remove maps that are too highly correlated
    Snets = zeros(params.N,params.N,params.S);
    for s = 1:params.S
        r = corr(P{s}); r = triu(r,1);
        Snets(:,:,s) = r;
    end
    Snets = mean(Snets,3);
    x = find(abs(Snets)>0.2); [x,y] = ind2sub(size(Snets),x);
    Modes_remove = unique(x);
    if length(Modes_remove)<10;
        X = setdiff(1:params.N,Modes_remove);
        Nextra = 10-length(Modes_remove);
        N = randperm(length(X)); N = N(1:Nextra);
        Modes_remove = [Modes_remove; X(N)'];
    end
    Modes_keep = setdiff(1:params.N,Modes_remove);
    for s = 1:params.S
        P{s} = P{s}(:,Modes_keep);
    end
    params.N = 15;
    
    if length(Modes_keep) == params.N
        
        fprintf('Succeeded (%d)\n',n);
        
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
        
        %     if plot_SamFigures
        %         cP = P{1}' * P{1};
        %         cP = cP ./ sqrt(diag(cP)*diag(cP)');
        %         figure; imagesc(cP, [-1 1]); colorbar; colormap(bluewhitered)
        %         xlabel('Mode'); ylabel('Mode');
        %         axis square
        %         set(gcf, 'Position', [200 200 500 400])
        %
        %         cA = [An{1}{1} An{1}{2} An{1}{3} An{1}{4}] * [An{1}{1} An{1}{2} An{1}{3} An{1}{4}]';
        %         cA = cA ./ sqrt(diag(cA)*diag(cA)');
        %         figure; imagesc(cA, [-1 1]); colorbar; colormap(bluewhitered)
        %         xlabel('Mode'); ylabel('Mode');
        %         axis square
        %         set(gcf, 'Position', [200 200 500 400])
        %     end
        
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
        
        %% Save nifti's and run latest versions of melodic and profumo
        
        Save_niftis(D,params,fileNameN,atlasParams);
        system(sprintf('sh Overlap_functions/ICA_PROFUMO.sh %1.2f %s %d %s',params.TR,fileNameN,params.iN,pwd))
        
        %     %% Extract PFMs
        %
        %     %For the first run, may want to plot convergence
        %     if plotNow
        %         %If we want to run with plots of convergence, pass in best guesses
        %         %for mean map and node time courses
        %         [pfmP1, pfmA1, pfmPg1] = runVBGP(D, params, params.iN, Pg, A);
        %     else
        %         %Just run VBGP normally
        %         [pfmP1, pfmA1, pfmPg1] = runVBGP(D, params, params.iN);
        %     end
        %
        %
        %     %Save scores
        %     scores.PFMs.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        %         makePA(pfmP1,pfmA1,params), params);
        %     [scores.PFMs.P(:,2*(n-1)+1), scores.PFMs.A(:,2*(n-1)+1), ...
        %         scores.PFMs.cP(:,2*(n-1)+1), scores.PFMs.cA(:,2*(n-1)+1)] ...
        %         = calculateDecompositionAccuracy(P, pfmP1, A, pfmA1, params);
        %
        %
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
        
        %Traditional dual regression
        addpath /usr/local/fmrib/fmt/
        tica_P1_DRnew = cell(params.S,1);
        tica_A1_DRnew = cell(params.S,1);
        for s = 1:params.S
            tica_A1_DRnew{s} = cell(params.R(s),1);
            Ds = zeros(params.V, params.T*params.R(s));
            M = zeros(atlasParams.V,params.iN,params.R(1));
            for r = 1:params.R(s)
                Ds(:, (r-1)*params.T+(1:params.T)) = D{s}{r};
                tica_A1_DRnew{s}{r} = demean((pinv(demean(double(ticaP1{1})))*demean(double(D{s}{r})))');
                M(:,:,r) = demean((pinv(demean(double(tica_A1_DRnew{s}{r})))*demean(double(D{s}{r}))')');
            end
            tica_P1_DRnew{s} = mean(M,3);
        end
        
        %Save scores
        scores.tICA_DR.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
            makePA(ticaP1_DR,ticaA1_DR,params), params);
        [scores.tICA_DR.P(:,2*(n-1)+1), scores.tICA_DR.A(:,2*(n-1)+1), ...
            scores.tICA_DR.cP(:,2*(n-1)+1), scores.tICA_DR.cA(:,2*(n-1)+1)] ...
            = calculateDecompositionAccuracy(P, ticaP1_DR, A, ticaA1_DR, params);
        
        
        %% Run Steve's spatio-temporal mixed ICA
        addpath ~steve/matlab/
        addpath ~steve/NETWORKS/FSLNets;
        
        % Concatenate data
        Dmat = zeros(params.V, params.T*sum(params.R));
        sr = 1;
        for s = 1:params.S
            for r = 1:params.R(s)
                Dmat(:, (sr-1)*params.T+(1:params.T)) = D{s}{r};
                sr = sr + 1;
            end
        end
        
        % Run spatio-temporal mixed ICA
        [pcaU,pcaS,pcaV]=nets_svds(Dmat,params.iN);
        Y1=pcaU'; Y2=pcaV';
        [A1,S1,A2,S2,A12,S12] = mfastica(Y1,Y2,10);
        
        % Rearrange timecourses
        mica_P1_DRnew = cell(params.S,1);
        micaA1 = cell(params.S,1); sr = 1;
        for s = 1:params.S
            micaA1{s} = cell(params.R(s),1);
            M = zeros(atlasParams.V,params.iN,params.R(1));
            for r = 1:params.R(s)
                micaA1{s}{r} = S12(:, (sr-1)*params.T+(1:params.T)+10000);
                M(:,:,r) = demean((pinv(demean(double(micaA1{s}{r})'))*demean(double(D{s}{r}))')');
                sr = sr + 1;
            end
            mica_P1_DRnew{s} = mean(M,3);
        end
        clear Dmat sr s r M
        
        % Get group maps
        micaP1 = cell(1,1);
        micaP1{1} = (A12'*Y1)';
        
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
        
        %% Save results
        clear Ds pcaA pcaP pcaS pcaU pcaV S1 S2 S12 svdS svdU svdV Y1 Y2 A1 A12 A2
        save(sprintf('Results/PFMsims_atlas_%s',fileNameN),'-v7.3')
        n = n+1;
        
        %% Load results from external runs of melodic and profumo
        %     [sicaPg_new,sica_P1_DR_new,sica_A1_DR_new,sica_P1_DR_thres,sica_A1_DR_thres] = Melodic_DR(fileName,D,atlasParams,params);
        %     [pfmPg1_new,pfmP1_new,pfmA1_new] = loadNewPROFUMO(fileName,params);
        
        %% Plot results
        if plot_JanineFigures
            RunOverlapTests_plots
        end
        
    else
        fprintf('Failed\n');
        
    end
end

%% Plot results

if plot_SamFigures
    plotScores(scores, params);
    input('Press return to continue')
    close all
end

% When finished - run next script to produce figure for paper
Paper_figure


