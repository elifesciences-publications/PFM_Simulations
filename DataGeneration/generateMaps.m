function [ P, plotHandle ] = generateMaps(params, options, plotFigures)
%Generates spatial maps for fMRI simulations
%   Can generate a variety of spatial structures

if nargin == 2
    plotFigures = false;
end

%Decide what function to call based on what type of map we want
switch options.P.form
    
    case 'Additive'
        P = generateAdditiveMaps(params, options, plotFigures);
        
    otherwise
        error('Not a recognised form for P')
        
end

%Then add registration
switch options.P.registration.form
    
    case 'RandomSmooth'
        P = registerMaps_RandomSmooth(P, params, options, plotFigures);
        
    case 'Null'
        
    otherwise
        error('Not a recognised form for P')
        
end

%Finally, add a global rescaling such that maps are overall unit variance
vP = 0;
for s = 1:params.S
    vP = vP + var(P{s}(:));
end
vP = vP / params.S;
for s = 1:params.S
    P{s} = P{s} / sqrt(vP);
end

if plotFigures
    plotMaps(P, params, options);
end

plotHandle = @plotMaps;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ P ] = generateAdditiveMaps(params, options, plotFigures)
%Generate a set of maps where the subject differences are simply added to
%the global maps

%Generate global maps
switch options.Pg.form
    case 'SS'
        Pg = generateAdditiveMaps_PgSS(params, options, plotFigures);
    case 'Freq'
        Pg = generateAdditiveMaps_PgFreq(params, options, plotFigures);
    case 'BlockAtlas'
        Pg = generateAdditiveMaps_PgBlockAtlas(params, options, plotFigures);
    case 'RandBoxcar'
        Pg = generateAdditiveMaps_RandBoxcar(params, options, plotFigures);
    case 'BiasedBoxcar'
        Pg = generateAdditiveMaps_BiasedBoxcar(params, options, plotFigures);
    otherwise
        error('Not a recognised form for Pg')
end
%Normalise
Pg = Pg / std(Pg(:));

%Generate subject maps
switch options.Ps.form
    case 'Gaussian'
        Ps = generateAdditiveMaps_PsGaussian(params, options, Pg, plotFigures);
    case 'SS'
        Ps = generateAdditiveMaps_PsSS(params, options, Pg, plotFigures);
    case 'Null'
        Ps = generateAdditiveMaps_PsNull(params, options, Pg, plotFigures);
    otherwise
        error('Not a recognised form for Ps')
end


%Combine global and subject maps, rescaling as necessary
P = cell(params.S, 1);
for s = 1:params.S
    if std(Ps{s}(:)) ~= 0
        P{s} = Pg + ( options.P.PsPg * Ps{s} / std(Ps{s}(:)) );
    else
        P{s} = Pg;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Pg ] = generateAdditiveMaps_PgSS(params, options, plotFigures)
%Generates a set of global maps following a spike slab distribution
%
% options.Pg.rot - amount of random rotations to induce correlations
% options.Pg.p - proportion of non zero elements in map
% options.Pg.sigma - std of non zero elements
% options.Pg.epsilon - std of additive noise

%Random rotation matrix
Pg_Rot = eye(params.N) + options.Pg.rot*( rand(params.N) - 0.5 );

%Generate a set of random maps
Pg = options.Pg.sigma * randn(params.V, params.N);
%Rotate to induce correlations
Pg = Pg * Pg_Rot;
%Sparsify according to spike-slab
Pg = Pg .* (rand(params.V, params.N) < options.Pg.p);
%Add some noise
Pg = Pg  + options.Pg.epsilon * randn(params.V, params.N);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Pg ] = generateAdditiveMaps_PgFreq(params, options, plotFigures)
%Generates a set of sparse global maps with enhanced low frequency content
%
% options.Pg.rot - amount of random rotations to induce correlations
% options.Pg.fc - frequency cut off (Nyquist frequency normalised to 1)
% options.Pg.fAmp - amount to amplify frequencies less than fc by
% options.Pg.p - proportion of map entries to retain
% options.Pg.epsilon - std of additive noise

%Spatial frequencies - Nyquist set to 1
f = linspace(0, 1, params.V/2+1);
%Randomly generate frequency amplitudes
fPg = abs(randn(length(f), params.N));
fPg(f < options.Pg.fc, :) = options.Pg.fAmp * fPg(f < options.Pg.fc, :);
%fPg = abs([ groundTruth.Pg.fAmp*randn(N, sum(f<groundTruth.Pg.fc)) ...
%    randn(N, sum(f>=groundTruth.Pg.fc)) ]);

%Add random phase
fPg = fPg .* exp( 2*pi*rand(size(fPg))*sqrt(-1) );
%Add conjugate 'high f' content for FFT
fPg = [fPg; conj(fPg(end-2:-1:1, :))];
%Now invert
Pg = ifft( fPg );
Pg = abs(Pg) .* sign(angle(Pg));
%Random rotation matrix
Pg_Rot = eye(params.N) + options.Pg.rot*( rand(params.N) - 0.5 );
%Rotate to induce correlations
Pg = Pg * Pg_Rot;
%Assume maps are Gaussian and find expected value of percentile given by p
Pg = Pg ./ std(Pg(:));
mu = norminv([options.Pg.p/2 1-options.Pg.p/2], 0, std(Pg(:)));
%Sparsify
Pg( (Pg > mu(1)) & (Pg < mu(2)) ) = 0;
%Add some noise
Pg = Pg  + options.Pg.epsilon * randn(params.V, params.N);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Pg ] = generateAdditiveMaps_PgBlockAtlas(params, options, plotFigures)
%Generates an atlas where activations are grouped in blocks
%  Just divides the set of voxels into params.N non-overlapping regions.
%  Region sizes are drawn from a Dirichlet distribution
%
% options.Pg.widthPrecision - parameter controlling variability in width


%Sample the block lengths
% - Dirichlet distribution
%This is a normalised vector of unit-scale gamma random variables
lengths = gamrnd(options.Pg.widthPrecision, 1, params.N, 1);
lengths = lengths / sum(lengths);
%Convert to voxels, taking into account sparsity
lengths = params.V * lengths;

%If we wanted to have the subject atlases varying around these lengths can
%use the following to generate subject lengths around the length mean
%lengthsS = gamrnd(repmat(lengths, 1, params.S)*10000, 1, params.N, params.S);
%lengthsS = bsxfun(@times, lengthsS, 1./sum(lengthsS));

Pg = zeros(params.V, params.N);
for n = 1:params.N
    inds = (1 + round(sum( lengths(1:(n-1)) ))):round(sum( lengths(1:n) ));
    Pg(inds, n) = 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Pg ] = generateAdditiveMaps_BiasedBoxcar(params, options, plotFigures)
%Generates a set of spatial maps where activations are grouped in blocks
%  A random number of blocks and a random sparsity map sparsity are
%  simulated. Then, the blocks lengths, signs and positions are generated.
%  Positions are set up such that there is a slight bias against overlap.
%  This is repeated to give a set of block based activations for each
%  network. The weights are drawn from a gamma distribution. Finally, the 
%  maps are smoothed slightly.
%
% options.P.nBlocks - expected number of separate blocks per map
% options.P.p - expected sparsity of the maps
% options.P.pVar - variance of this sparsity over networks
% options.P.pPosBlock - proportion of blocks that are positive
% options.P.smootherWidth - width, in voxels, of the smoothing filter
% options.P.minWeight - allows a minimum map weight to be specified
% options.P.weightRange.a - gamma weight shape parameter
% options.P.weightRange.b - gamma weight rate parameter
% options.P.biasStrength - how strong the exclusion bias is, in [0 1]

%Random sparsity (overall number of voxels in the blocks) is Beta distributed
%First recover the a and b parameters from the mean and variance
m = options.P.p; v = options.P.pVar;
a = (m^2 * (1-m) / v) - m;
b = (m * (1-m)^2 / v) - (1-m);

if plotFigures
    %Plot the sparsity distribution
    x = linspace(0, 1, 250);
    p = betapdf(x, a, b);
    figure; plot(x, p);
    xlim([0 1]); xlabel('x'); ylabel('p(x)')
    title('BiasedBoxcar: distribution of map sparsity')
    
    %Plot the weight distribution
    gamMode = max(1, (options.P.weightRange.a-1)/options.P.weightRange.b);
    x = linspace(0, 5*gamMode, 500);
    p = gampdf(x, options.P.weightRange.a, 1/options.P.weightRange.b);
    
    figure; plot(x+options.P.minWeight, options.P.p * options.P.pPosBlock * p);
    hold on; plot(-x-options.P.minWeight, options.P.p * (1-options.P.pPosBlock) * p);
    plot([0 0], [0 (1-options.P.p)], 'r');
    xlim([-x(end) x(end)]); xlabel('x'); ylabel('p(x)')
    title('BiasedBoxcar: distribution of map weights')
    
end

%Generate empty maps
Pblocks = zeros(params.V, params.N);

%Generate the parameters of the global distribution
for n = 1:params.N
    
    %Generate a random number of blocks
    % - Poisson distributed)
    nBlocks = poissrnd(options.P.nBlocks) + 1;
    
    %Generate a random sparsity parameter from the distribution
    p = betarnd(a,b);
    
    %Sample the block lengths
    % - Dirichlet distribution
    %This is a normalised vector of unit-scale gamma random variables
    %First make the dirichlet prior - by making the first element half the
    %total this controls the minimum block size
    lengths = ones(nBlocks,1); lengths(1) = max(1, 1.0*(nBlocks-1));
    lengths = lengths(randperm(nBlocks));
    lengths = gamrnd(10*lengths, 1, nBlocks, 1); %First param controls variability
    lengths = lengths / sum(lengths);
    %Convert to voxels, taking into account sparsity
    lengths = p * params.V * lengths;
    lengths = round(lengths); lengths(lengths==0) = 1;
    
    %Sample the block start points
    if n == 1
        %Sample first set randomly
        startPoints = randi(params.V - max(lengths) + 1, nBlocks);
    else
        %For subsequent blocks, bias against spatial overlap
        %Find out how overlapping previous blocks are
        spatialOverlap = sum(abs(Pblocks(:, 1:(n-1))), 2);
        spatialOverlap = options.P.biasStrength * (spatialOverlap) ...
            + (1 - options.P.biasStrength) * mean(spatialOverlap);
        
        startPoints = zeros(nBlocks,1);
        for b1 = 1:nBlocks
            %Generate a weighting that biases against starting somewhere
            %where overlap is high
            weights = conv(1./spatialOverlap, ones(lengths(b1),1), 'valid');
            if sum(isinf(weights) ~= 0)
                weights(~isinf(weights)) = 0; weights(isinf(weights)) = 1;
            end
            inds = 1:(params.V-lengths(b1)+1);
            %Randomly sample a start point
            startPoints(b1) = randsample(inds, 1, true, weights(inds));
        end
    end
    
    %Sample the block signs
    blockSigns = sign(options.P.pPosBlock - rand(nBlocks, 1));
    
    %Sample the block weights
    %blockWeights = options.P.minWeight ...
    %    + gamrnd(options.P.weightRange.a, 1/options.P.weightRange.b, nBlocks, 1);
    
    %Add the blocks in
    for bl = 1:nBlocks
        endPoint = min(params.V, startPoints(bl)+lengths(bl)-1);
        %Pblocks(startPoints(bl):endPoint, n) = blockSigns(bl) * blockWeights(bl);
        Pblocks(startPoints(bl):endPoint, n) = blockSigns(bl) ...
            * (gamrnd(options.P.weightRange.a, 1/options.P.weightRange.b, ...
            lengths(bl), 1) + options.P.minWeight);
    end
    
end

%Finally, smooth the block maps with a box filter
h = ones(options.P.smootherWidth,1) / options.P.smootherWidth;
for n = 1:params.N
    Pblocks(:,n) = conv(Pblocks(:,n), h, 'same');
end

Pg = Pblocks;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Ps ] = generateAdditiveMaps_PsGaussian(params, options, Pg, plotFigures)
%Generates subject maps following a Gaussian distribution

Ps = cell(params.S, 1);
for s = 1:params.S
    Ps{s} = randn(params.V, params.N);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Ps ] = generateAdditiveMaps_PsSS(params, options, Pg, plotFigures)
%Generates subject maps following a spike slab distribution (in reality
%repeatedly calls generateAdditiveMaps_PgSS)
%
% options.Ps.rot - amount of random rotations to induce correlations
% options.Ps.p - proportion of non zero elements in map
% options.Ps.sigma - std of non zero elements
% options.Ps.epsilon - std of additive noise

%Rename options for global function
PsOptions.Pg.rot = options.Ps.rot;
PsOptions.Pg.p = options.Ps.p;
PsOptions.Pg.sigma = options.Ps.sigma;
PsOptions.Pg.epsilon = options.Ps.epsilon;

%Call global function with new options to get subject maps
Ps = cell(params.S, 1);
for s = 1:params.S
    Ps{s} = generateAdditiveMaps_PgSS(params, PsOptions, plotFigures);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Ps ] = generateAdditiveMaps_PsNull(params, options, Pg, plotFigures)
%Generates a null set of subject maps i.e. all zero

Ps = cell(params.S, 1);
for s = 1:params.S
    Ps{s} = zeros(params.V, params.N);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ P ] = registerMaps_RandomSmooth(P, params, options, plotFigures)
%Modifies input spatial maps by a smooth warp
% Smooth means voxels cannot chnage order in the warp
%
% options.P.registration.maxError - max possible deviation, in voxels

%Loop over all subjects, generating new registrations each time
for s = 1:params.S
    
    %Generate a random vector and convolve with a boxcar to make smooth
    % This becomes the rate of change of the warp
    R = randn(1, params.V); L = floor(0.05*params.V);
    Rs = conv(R, ones(L,1)/L, 'same');
    Rs = Rs - mean(Rs);
    
    if plotFigures && (s==1)
        figure; plot(R, '.');% hold on; plot(Rs, 'r')
    end
    
    %Normalise size of Rs and squash into -1 < Rs < 1
    Rs = Rs / std(Rs);
    Rs(Rs<0) = -( 1 - exp(Rs(Rs<0)) );
    Rs(Rs>0) = 1 - exp(-Rs(Rs>0));
    Rs = Rs - mean(Rs);
    
    if plotFigures && (s==1)
        %Plot the warp gradient
        hold on; plot(Rs, 'r')
        xlim([1 params.V]); xlabel('Voxel')
        legend('Random noise', 'Smoothed warp gradient')
        title('Registration: random warp gradient')
    end
    
    %Now integrate this rate of change to get the actual warp
    %Conditions on Rs are:
    % - Zero mean (no overall drift)
    % - Rs > -1 (no change in order)
    W = cumtrapz(Rs);
    W = options.P.registration.maxError * W / max(abs(W));
    
    if plotFigures && (s==1)
        %Plot the warp
        figure; plot(1:params.V, 1:params.V, '--');
        hold on; plot(1:params.V, (1:params.V)+W, 'r');
        xlim([1 params.V]); xlabel('Voxel in')
        ylim([1 params.V]); ylabel('Voxel out')
        title('Registration: random warp')
    end
    
    %Use the warp interpolate the spatial maps
    for n = 1:params.N
        P{s}(:,n) = interp1(1:params.V, P{s}(:,n), (1:params.V) + W, 'pchip');
    end
    
    %if plotFigures && (s==1)
    %    %Plot an example of the interpolation
    %    figure; plot((1:0.1:100), interp1(1:params.V, P{s}(:,n), (1:0.1:100), 'pchip', 0))
    %    hold on; plot(1:100, P{s}(1:100,n), 'r+')
    %end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = plotMaps(P, params, options)

%Plot the mean map
Pg = 0;
for s = 1:params.S
    Pg = Pg + P{s};
end
Pg = Pg / params.S;
figure; imagesc(Pg, max(abs(Pg(:)))*[-1 1]);
colorbar; colormap(bluewhitered)
title('Mean map')

%And the distribution of weights
figure; TwoColourTufteHist(Pg(:), 'normalise')
title(sprintf('Mean map distribution (kurtosis: %.2f)', kurtosis(Pg(:))))

%Spatial correlations
cP = Pg'*Pg;
dcP = sqrt( diag(cP) );
cP = cP ./ (dcP * dcP');
figure; imagesc(cP, [-1 1]); colorbar; colormap(bluewhitered)
title('Correlations in mean map')

%Amount of overlap
figure; plot(sum(abs(Pg),2));
xlim([1 params.V]); xlabel('Voxel')
ylabel('Sum, over networks, of mean map')
title('Fluctuations in spatial overlap')

%Plot one subject map
s = randi(params.S,1);
figure; imagesc(P{s}, max(abs(Pg(:)))*[-1 1]);
colorbar; colormap(bluewhitered)
title(['Subject ' num2str(s) ' Map'])

%Plot all the subject versions of one network
Ps = NaN(params.V, params.S+1); n = randi(params.N,1);
for s = 1:params.S
    Ps(:,s) = P{s}(:,n);
end
Ps(:,end) = Pg(:,n);
figure; imagesc(Ps, max(abs(Pg(:)))*[-1 1]);
colorbar; colormap(bluewhitered)
%export_fig('SubjectMaps', '-pdf', '-nocrop', '-transparent')
title(['Subject versions of network ' num2str(n) ' map'])

% Look at how similar subject maps are to each other and the mean map
cP = Pg;
for s = 1:params.S
    cP = [cP P{s}];
end
cP = cP' * cP;
dcP = sqrt( diag(cP) );
cP = cP ./ (dcP * dcP');

%First to the mean map
PsPg = NaN(params.N, params.S); inds = 1:params.N;
for s = 1:params.S
    PsPg(:,s) = diag(cP( inds, s*params.N+inds ));
end
figure; TwoColourTufteHist(PsPg(:), 'xlim', [-1 1], 'normalise');
xlabel('Correlation'); ylabel('Relative Frequency')
title('Subject map - mean map correlations')

%Then to each other
cP = cP((params.N+1):end, (params.N+1):end); PsPs = [ ];
for s = 1:(params.S-1)
    PsPs = [PsPs; diag(cP, s*params.N)];
end
figure; TwoColourTufteHist(PsPs(:), 'xlim', [-1 1], 'normalise');
xlabel('Correlation'); ylabel('Relative Frequency')
title('Subject map - subject map correlations')

input('Press return to continue')
close all

end
