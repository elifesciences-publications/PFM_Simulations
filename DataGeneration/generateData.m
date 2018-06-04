function [ D, plotHandle ] = generateData(PA, params, options, plotFigures)
%Generates noisy scans given observed BOLD signal

if nargin == 3
    plotFigures = false;
end

%Decide what function to call based on what type of timecourse we want
switch options.D.noise.form
    
    case 'SpatiotemporallyWhiteGaussian'
        D = generateData_STWhiteGaussian(PA, params, options, plotFigures);
        
    case 'SpatiotemporallyWhiteTdist'
        D = generateData_STWhiteTdist(PA, params, options, plotFigures);
        
    case 'StructuredTdist'
        D = generateData_StructuredTdist(PA, params, options, plotFigures);
        
    otherwise
        error('Not a recognised form for D')
        
end

if plotFigures
    plotData(D, params, options);
end

plotHandle = @plotData;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ D ] = generateData_STWhiteGaussian(PA, params, options, plotFigures)
%Adds independent Gaussian noise to everything to match required SNR

%Find global variance of BOLD signal
vPA = 0;
for s = 1:params.S
    for r = 1:params.R(s)
        vPA = vPA + var(PA{s}{r}(:));
    end
end
vPA = vPA / sum(params.R);

%Use SNR to find noise std from signal var
% SNR = var(signal) / var(noise)
noiseStd = sqrt( vPA / options.D.SNR );

%Add Gaussian noise to the signal
D = cell(params.S, 1);
for s = 1:params.S
    D{s} = cell(params.R(s), 1);
    for r = 1:params.R(s)
        D{s}{r} = PA{s}{r} + noiseStd * randn(params.V, params.T);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ D ] = generateData_STWhiteTdist(PA, params, options, plotFigures)
%Adds independent t-distributed noise to everything to match required SNR
%
% options.D.noise.a - gamma precision shape parameter
% options.D.noise.b - gamma precision rate parameter

%Plot probability distribution of precision parameters
if plotFigures
    gamMode = max(1, (options.D.noise.a-1)/options.D.noise.b);
    x = linspace(0, 3*gamMode, 500);
    
    figure; plot(x, gampdf(x, options.D.noise.a, 1/options.D.noise.b));
    xlabel('x'); xlim([x(1) x(end)]); ylabel('p(x)')
    title('t-distribution precisions')
end

%Find global variance of BOLD signal
vPA = 0;
for s = 1:params.S
    for r = 1:params.R(s)
        vPA = vPA + var(PA{s}{r}(:));
    end
end
vPA = vPA / sum(params.R);

%Use SNR to find noise std from signal var
% SNR = var(signal) / var(noise)
noiseStd = sqrt( vPA / options.D.SNR );

%Add Gaussian noise to the signal
D = cell(params.S, 1);
for s = 1:params.S
    D{s} = cell(params.R(s), 1);
    for r = 1:params.R(s)
        precisions = gamrnd(options.D.noise.a, 1/options.D.noise.b, params.V, params.T);
        noise = randn(params.V, params.T) ./ sqrt(precisions);
        
        D{s}{r} = PA{s}{r} + noiseStd * noise / std(noise(:));
    end
end

%Plot noise distribution
if plotFigures
    figure; TwoColourTufteHist(noiseStd * noise / std(noise(:)), 'normalise')
    title(sprintf('Noise distribution (kurtosis: %.2f)', kurtosis(noise(:))))
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ D ] = generateData_StructuredTdist(PA, params, options, plotFigures)
%Adds:
%  + a low rank structured subspace
%  + independent t-distributed noise
%
% options.D.noise.a - gamma precision shape parameter
% options.D.noise.b - gamma precision rate parameter
% options.D.noise.structuredStd - std of structured subspace
% options.D.noise.unstructuredStd - std of unstructured subspace
% options.D.noise.N - rank of structured noise

%Find global variance of BOLD signal
vPA = 0;
for s = 1:params.S
    for r = 1:params.R(s)
        vPA = vPA + var(PA{s}{r}(:));
    end
end
vPA = vPA / sum(params.R);

%Use SNR to find noise std from signal var
% SNR = var(signal) / var(noise)
noiseStd = sqrt( vPA / options.D.SNR );

%Add noise to the signal
D = cell(params.S, 1);
for s = 1:params.S
    D{s} = cell(params.R(s), 1);
    for r = 1:params.R(s)
        % Low rank subspace
        structured = randn(params.V, options.D.noise.N) * randn(options.D.noise.N, params.T);
        structured = options.D.noise.structuredStd * structured / std(structured(:));
        
        % T-distribution, unstructured
        precisions = gamrnd(options.D.noise.a, 1/options.D.noise.b, params.V, params.T);
        unstructured = randn(params.V, params.T) ./ sqrt(precisions);
        unstructured = options.D.noise.unstructuredStd * unstructured / std(unstructured(:));
        
        noise = structured + unstructured;
        D{s}{r} = PA{s}{r} + noiseStd * noise / std(noise(:));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = plotData(D, params, options)

%Plot a couple of example scans
s = randi(params.S,1); r = randi(params.R(s),1);
figure; imagesc(D{s}{r}, 6*[-1 1]); colorbar
title(['Example scan: subject ' num2str(s) ', repeat ' num2str(r)])

s = randi(params.S,1); r = randi(params.R(s),1);
figure; imagesc(D{s}{r}, prctile(abs(D{s}{r}(:)),100)*[-1 1]); colorbar
title(['Example scan: subject ' num2str(s) ', repeat ' num2str(r)])

input('Press return to continue')
close all

end
