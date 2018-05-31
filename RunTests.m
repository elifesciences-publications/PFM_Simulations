%Generates a test fMRI data set and compares the performance of a variety
%of node discovery algorithms

% fsl_sub -q verylong.q matlab -singleCompThread -nojvm -nosplash -nodisplay \< RunTests.m
% samh@jalapeno18 $ /usr/local/MATLAB/R2012a/bin/matlab -nojvm -nosplash -nodisplay -r RunTests

restoredefaultpath
addpath(genpath('~/Documents/Code/MATLAB/'))
addpath(genpath('~/Documents/Code/Algorithms/'))
addpath('~/Documents/Code/VBGP_Tests/')
addpath('DataGeneration/'); addpath('Methods/'); 
addpath('Scoring/'); addpath('Visualisation/')

close all; clear all; clc

rng('shuffle')

fileName = '~/scratch/SimulatedTests'

plotFigures = false;
if plotFigures
    prettyFigures()
end

%% Set size of problem

%Number of times to repeat simulation / test cycle
params.nRepeats = 10;

%Details of scans
params.S = 30;       %Subjects
params.R = 4*ones(params.S,1);   %Repeats

params.T = 1200;     %No. of time points per fMRI scan
params.TR = 0.72;
params.dt = 0.1;     %Neural sampling rate
%Amount of neural points to simulate - more than scan length so don't have
%to zero pad HRF at start of scan
params.Tn = ceil(1.25 * params.T * params.TR / params.dt);

params.V = 7500;     %Voxels

params.N = 25;       %Number of nodes


%% Set the details of the tests we want to do

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spatial maps
options.P.form = 'Additive';

%Choose form for Pg
options.Pg.form = 'RandBoxcar';
%Set options based on that choice
switch options.Pg.form
    
    case 'SS'
        options.Pg.rot = 1;
        options.Pg.p = 0.2;
        options.Pg.sigma = 1;
        options.Pg.epsilon = 0.1;
        
    case 'Freq';
        options.Pg.rot = 0.5;
        options.Pg.p = 0.2;
        options.Pg.fc = 0.01;
        options.Pg.fAmp = 10;
        options.Pg.epsilon = 0.1;
        
    case 'BlockAtlas';
        options.Pg.widthPrecision = 100;
        
    case 'RandBoxcar';
        options.P.nBlocks = 4;
        options.P.p = 0.1;
        options.P.pVar = 0.0005;
        options.P.pPosBlock = 0.5;
        options.P.smootherWidth = floor( 0.01*params.V );
        options.P.weightRange = 0.5;
end

%Choose form for Ps
options.Ps.form = 'SS';
options.P.PsPg = 0.75; %Ratio of std(P{s}(:)) / std(Pg(:))
%Set options based on that choice
switch options.Ps.form
    
    case 'SS'
        options.Ps.rot = 1;
        options.Ps.p = 0.1;
        options.Ps.sigma = 1;
        options.Ps.epsilon = 0.1;
        
    case 'Gaussian'
        
    case 'Null'
end

%Choose registration errors
options.P.registration.form = 'RandomSmooth';
%options.P.registration.maxError = 0.3 * params.V / params.N;
options.P.registration.maxError = 0.015 * params.V;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Time courses
options.An.form = 'Freq';
switch options.An.form
    case 'Freq';
        options.An.rot = 0.75;
        options.An.rotS = 0.25;
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
        %options.BS.satCoeff = 0.075;
        options.BS.tanhPercentile = 99;
        options.BS.tanhMax = 1.25;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Noise
options.D.noise.form = 'SpatiotemporallyWhiteGaussian';
%Signal to noise ratio (expressed in terms of power)
options.D.SNR = 0.1;

%% Run the tests

for n = 1:params.nRepeats
    
    %% Generate data
    if plotFigures && (n==1)
        plotNow = true;
    else
        plotNow = false;
    end
    
    P = generateMaps(params, options, plotNow);
    
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
    
    [scores.GTg.P(:,n), scores.GTg.A(:,n)] = calculateDecompositionAccuracy( ...
        P, repmat({Pg}, params.S, 1), A, Ag, params);
    
    %% Run VBGP
    
    %For the first run, may want to plot convergence
    if plotNow
        %If we want to run with plots of convergence, pass in best guesses
        %for mean map and node time courses
        [vbgpP1, vbgpA1] = runVBGP(D, params, params.N, Pg, A);
    else
        %Just run VBGP normally
        [vbgpP1, vbgpA1] = runVBGP(D, params, params.N);
    end
    
    
    %Save scores
    scores.VBGP.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(vbgpP1,vbgpA1,params), params);
    [scores.VBGP.P(:,2*(n-1)+1), scores.VBGP.A(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, vbgpP1, A, vbgpA1, params);
    
    
    
    %Repeat analysis
    [vbgpP2, vbgpA2] = runVBGP(D, params, params.N);
    
    %Save scores
    scores.VBGP.PA(2*(n-1)+2) = calculateBoldRecovery(PA, ...
        makePA(vbgpP2,vbgpA2,params), params);
    [scores.VBGP.P(:,2*(n-1)+2), scores.VBGP.A(:,2*(n-1)+2)] ...
        = calculateDecompositionAccuracy(P, vbgpP2, A, vbgpA2, params);
    
    
    
    % Test - retest
    [scores.VBGP.Ptr(:,2*(n-1)+1), scores.VBGP.Ptr(:,2*(n-1)+2)] ...
        = calculateTestRetestScores(vbgpP1, vbgpP2, params);
    
    %% Run PCA (using SVD)
    
    [pcaP, pcaA, svdU, svdS, svdV] = runPCA(D, params, params.N);
    
    %fprintf('Hello, Sam!! \n');
    
    %Save scores
    scores.PCA.PA(n) = calculateBoldRecovery(PA, ...
        makePA(repmat({pcaP},params.S,1),pcaA,params), params);
    
    [scores.PCA.P(:,n), scores.PCA.A(:,n)] = calculateDecompositionAccuracy( ...
        P, repmat({pcaP}, params.S, 1), A, pcaA, params);
    
    %% Run sICA on low dim PCA, with dual reg
    
    [sicaP1, sicaA1] = runSICA(svdU, svdS, svdV, params, params.N);
    sicaP1 = repmat({sicaP1},params.S,1);
    
    %Save scores
    scores.sICA.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(sicaP1,sicaA1,params), params);
    [scores.sICA.P(:,2*(n-1)+1), scores.sICA.A(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, sicaP1, A, sicaA1, params);
    
    %Run dual regression
    [sicaP1_DR, sicaA1_DR] = runDR(D, sicaA1, params);
    
    %Save scores
    scores.sICA_DR.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(sicaP1_DR,sicaA1_DR,params), params);
    [scores.sICA_DR.P(:,2*(n-1)+1), scores.sICA_DR.A(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, sicaP1_DR, A, sicaA1_DR, params);
    
    
    
    %Repeat
    [sicaP2, sicaA2] = runSICA(svdU, svdS, svdV, params, params.N);
    sicaP2 = repmat({sicaP2},params.S,1);
    
    %Save scores
    scores.sICA.PA(2*(n-1)+2) = calculateBoldRecovery(PA, ...
        makePA(sicaP2,sicaA2,params), params);
    [scores.sICA.P(:,2*(n-1)+2), scores.sICA.A(:,2*(n-1)+2)] ...
        = calculateDecompositionAccuracy(P, sicaP2, A, sicaA2, params);
    
    %Run dual regression
    [sicaP2_DR, sicaA2_DR] = runDR(D, sicaA2, params);
    
    %Save scores
    scores.sICA_DR.PA(2*(n-1)+2) = calculateBoldRecovery(PA, ...
        makePA(sicaP2_DR,sicaA2_DR,params), params);
    [scores.sICA_DR.P(:,2*(n-1)+2), scores.sICA_DR.A(:,2*(n-1)+2)] ...
        = calculateDecompositionAccuracy(P, sicaP2_DR, A, sicaA2_DR, params);
    
    
    
    % Test - retest
    [scores.sICA.Ptr(:,2*(n-1)+1), scores.sICA.Ptr(:,2*(n-1)+2)] ...
        = calculateTestRetestScores(sicaP1, sicaP2, params);
    [scores.sICA_DR.Ptr(:,2*(n-1)+1), scores.sICA_DR.Ptr(:,2*(n-1)+2)] ...
        = calculateTestRetestScores(sicaP1_DR, sicaP2_DR, params);
    
    %% Run tICA on low dim PCA, with dual reg
    
    [ticaP1, ticaA1] = runTICA(svdU, svdS, svdV, params, params.N);
    ticaP1 = repmat({ticaP1},params.S,1);
    
    %Save scores
    scores.tICA.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(ticaP1,ticaA1,params), params);
    [scores.tICA.P(:,2*(n-1)+1), scores.tICA.A(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, ticaP1, A, ticaA1, params);
    
    %Run dual regression
    [ticaP1_DR, ticaA1_DR] = runDR(D, ticaA1, params);
    
    %Save scores
    scores.tICA_DR.PA(2*(n-1)+1) = calculateBoldRecovery(PA, ...
        makePA(ticaP1_DR,ticaA1_DR,params), params);
    [scores.tICA_DR.P(:,2*(n-1)+1), scores.tICA_DR.A(:,2*(n-1)+1)] ...
        = calculateDecompositionAccuracy(P, ticaP1_DR, A, ticaA1_DR, params);
    
    
    
    %Repeat
    [ticaP2, ticaA2] = runTICA(svdU, svdS, svdV, params, params.N);
    ticaP2 = repmat({ticaP2},params.S,1);
    
    %Save scores
    scores.tICA.PA(2*(n-1)+2) = calculateBoldRecovery(PA, ...
        makePA(ticaP2,ticaA2,params), params);
    [scores.tICA.P(:,2*(n-1)+2), scores.tICA.A(:,2*(n-1)+2)] ...
        = calculateDecompositionAccuracy(P, ticaP2, A, ticaA2, params);
    
    %Run dual regression
    [ticaP2_DR, ticaA2_DR] = runDR(D, ticaA2, params);
    
    %Save scores
    scores.tICA_DR.PA(2*(n-1)+2) = calculateBoldRecovery(PA, ...
        makePA(ticaP2_DR,ticaA2_DR,params), params);
    [scores.tICA_DR.P(:,2*(n-1)+2), scores.tICA_DR.A(:,2*(n-1)+2)] ...
        = calculateDecompositionAccuracy(P, ticaP2_DR, A, ticaA2_DR, params);
    
    
    
    % Test - retest
    [scores.tICA.Ptr(:,2*(n-1)+1), scores.tICA.Ptr(:,2*(n-1)+2)] ...
        = calculateTestRetestScores(ticaP1, ticaP2, params);
    [scores.tICA_DR.Ptr(:,2*(n-1)+1), scores.tICA_DR.Ptr(:,2*(n-1)+2)] ...
        = calculateTestRetestScores(ticaP1_DR, ticaP2_DR, params);
    
    %% Save scores to file
    
    save(fileName, 'scores', 'params', 'options')
    
    %%
    % P = randn(V,N);
    % figure; hold all;
    % for n = 1:50
    % n
    % score = paircomponents(P,Pg); plot(abs(score(:,3)), 'Color', [0 0 1] + 0.5*(1-n/50)*[1 1 0]);
    % pause(1)
    % A = inv(P'*P)*P'*Dmat;
    % P = Dmat*A'*inv(A*A');
    % end
    
    
end

%% Plot results

if plotFigures
    plotScores(scores, params);
    input('Press return to continue')
    close all
end

%% If run from the command line make sure we quit MATLAB

exit