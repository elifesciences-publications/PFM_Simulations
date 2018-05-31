function [  ] = plotScores( scores, params, saveFigures )
%Plots the results from whichever analyses have been run

if nargin == 2
    saveFigures = false;
end

if saveFigures
    figPosLong = 200 + [0 0 1300 200];
    figPosSquare = 200 + [0 0 300 300];
end
ticks = 0:0.2:1; %Where to place ticks for correlations

%% Plot the map decomposition accuracy

%Plot map accuracy (P)
[Pmethods, Presults] = extractResults(scores, params, 'P');

%As a boxplot
makeBoxPlot( Pmethods, Presults );
ylim([0 1]+0.05*[-1 1]); set(gca, 'YTick', ticks);
ylabel('Correlation coefficient')
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_SpatialMaps', '-pdf', '-transparent')
end
title('Spatial Map Accuracy')

%And all the values
plotAllVals(Pmethods, Presults, [0 1]+0.05*[-1 1])
set(gca, 'YTick', ticks)
ylabel('Correlation coefficient')
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_SpatialMapsAll', '-pdf', '-transparent')
end
title('Spatial Map Accuracy')

%% Plot the temporal accuracy

%Plot temporal results (A)
[Amethods, Aresults] = extractResults(scores, params, 'A');

makeBoxPlot( Amethods, Aresults );
ylim([0 1]+0.05*[-1 1]); set(gca, 'YTick', ticks);
ylabel('Correlation coefficient')
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_TimeCourses', '-pdf', '-transparent')
end
title('Time Course Accuracy')


%% Plot the spatial correlation accuracy

%Plot correlation results (cP)
[cPmethods, cPresults] = extractResults(scores, params, 'cP');

makeBoxPlot( cPmethods, cPresults );
ylabel('RMS Error'); ylims = get(gca, 'YLim'); ylims(1) = 0; ylim(ylims)
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_SpatialCorrs', '-pdf', '-transparent')
end
title('Spatial Correlation Accuracy')


%% Plot the temporal correlation accuracy

%Plot correlation results (cA)
[cAmethods, cAresults] = extractResults(scores, params, 'cA');

makeBoxPlot( cAmethods, cAresults );
ylabel('RMS Error'); ylims = get(gca, 'YLim'); ylims(1) = 0; ylim(ylims)
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_TemporalCorrs', '-pdf', '-transparent')
end
title('Temporal Correlation Accuracy')


%% Plot the ability to recover the latent space

%Plot ability to recover BOLD signal (PA)
[PAmethods, PAresults] = extractResults(scores, params, 'PA');

makeBoxPlot(PAmethods, PAresults);
ylabel('RMS Error'); ylims = get(gca, 'YLim'); ylims(1) = 0; ylim(ylims)
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_BOLD', '-pdf', '-transparent')
end
title('BOLD Signal Recovery')


%% Plot the test-retest scores

%Plot correlation results (Ptr)
[PTRmethods, PTRresults] = extractResults(scores, params, 'Ptr');

makeBoxPlot( PTRmethods, PTRresults );
ylim([0 1]+0.05*[-1 1]); set(gca, 'YTick', ticks);
ylabel('Correlation coefficient')
if saveFigures
    set(gcf, 'Position', figPosLong)
    export_fig('SimData_SpatialMapTR', '-pdf', '-transparent')
end
title('Spatial Map Repeatability')


%%
input('Press return to continue')


%% Plot the relationships between map and time course accuracy

range.x = [0 1]; range.y = [0 1];
H = makeScatterPlots(Pmethods, Presults, Amethods, Aresults, range);
for n = 1:length(H)
    set(0, 'CurrentFigure', H(n))
    axis square
    xlabel('Map accuracy')
    ylabel('Time course accuracy')
    set(gca, 'XTick', ticks); set(gca, 'YTick', ticks)
end

input('Press return to continue')
close(H)


%% Plot the relationships between spatial and temporal correlation accuracy

range.x = [0 0.7]; range.y = [0 0.7];
H = makeScatterPlots(cPmethods, cPresults, cAmethods, cAresults, range);
for n = 1:length(H)
    set(0, 'CurrentFigure', H(n))
    axis square
    xlabel('Spatial correlation accuracy')
    ylabel('Temporal correlation accuracy')
end

input('Press return to continue')
close(H)


%% Plot the relationship between map accuracy and test-retest

H = makeAccuracyRetestPlots(PTRmethods, PTRresults, Pmethods, Presults);
for n = 1:length(H)
    figure(H(n))
    set(0, 'CurrentFigure', H(n))
    set(gca, 'XTick', ticks); set(gca, 'YTick', ticks)
    if saveFigures
        method = get(gca,'Title'); method = get(method,'String'); title('')
        set(gcf, 'Position', figPosSquare)
        export_fig(['SimData_Acc-v-TR_' method], '-pdf', '-transparent')
        title(method)
    end
end

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ methods, results ] = extractResults( scores, params, testName )
%Extracts the set of scores for a given test from the 'scores' structure
%   Returns a list of the methods with that score, as well as the scores
%   themselves

%Extract the methods that have been tested
methods = fieldnames(scores);

%Loop over methods, extracting scores where appropriate
n = 1;
while n <= numel(methods)
    
    %See if this method has a PA score
    if isfield(scores.(methods{n}), testName)
        %If so, record it
        results{n} = scores.(methods{n}).(testName);
        %Move on to the next method
        n = n+1;
    else
        %If not present, remove the method
        methods = methods([1:(n-1) (n+1):end]);
    end
    
end

for n = 1:numel(methods)
    if strcmp(methods{n}, 'VBGP')
        methods{n} = 'PFMs';
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [  ] = makeBoxPlot( methods, results )
%Given a set of methods and results, generates a boxplot

%Find max number of data points recorded for any method
mLength = 0;
for n = 1:numel(results)
    mLength = max(mLength, numel(results{n}));
end

%Turn the results into a matrix for the boxplot function
boxmat = NaN(mLength, numel(methods));
for n = 1:numel(methods)
    boxmat(1:numel(results{n}), n) = results{n}(:);
end

%Plot results with appropriate labels
figure;
h = boxplot(boxmat, 'labels', methods);%, 'labelorientation', 'inline');
%Change line thicknesses
for ih=1:6
    set(h(ih,:),'LineWidth',2);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ H ] = makeScatterPlots( methods1, results1, methods2, results2, range )
%Draws scatter plots where the two sets of methods overlap

%Find common methods
[methods, inds1, inds2] = intersect(methods1, methods2);

%Loop over methods
H = [ ];
for n = 1:numel(methods)
    %If they have the same number of results do a scatter plot
    if numel(results1{inds1(n)}) == numel(results2{inds2(n)})
        
        h = figure; H = [H h]; hold on
        %If range has been specified, plot a line
        if nargin == 5
            %This just shows linear trend through the specified ranges
            plot(range.x, range.y, '--', 'Color', 0.8*[1 1 1])
        end
        %Plot the data
        plot(results1{inds1(n)}(:), results2{inds2(n)}(:), ...
            'b.', 'MarkerSize', 10);
        %Similarly set axes lims to the range
        if nargin == 5
            xlim(range.x); ylim(range.y);
        end
        title(methods{n}, 'interpreter', 'none')
        
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ] = plotAllVals( methods, results, lims )
%Plots all the results, rather than combining into a single box plot

methodSpacing = 1/2; repeatSpacing = 1/4;

figure; hold on; box on

xEnd = 0; labelTicks = NaN(numel(methods,1));
for n = 1:numel(methods)
    
    %Record where this method starts
    xStart = xEnd;
    x = xStart + methodSpacing;
    %Plot dividing line between methods
    if n ~= 1
        plot(xStart*[1 1], lims, '--', 'Color', 0.8*[1 1 1])
    end
    
    %Loop over methods, plotting all their points
    for r = 1:size(results{n},2);
        
        plot(x, results{n}(:,r), 'b.', 'MarkerSize', 10)
        x = x + repeatSpacing; %Advance to next result
        
    end
    
    %Record where we ended up, and where the label should go
    xEnd = (x - repeatSpacing) + methodSpacing;
    labelTicks(n) = (xStart + xEnd) / 2;
    
end

%Set lims
xlim([0 xEnd]); ylim(lims)

%Add labels
set(gca, 'XTick', labelTicks)
set(gca, 'XTickLabel', methods)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ H ] = makeAccuracyRetestPlots( PTRmethods, PTRresults, Pmethods, Presults )
%Draws scatter plots of the ground truth accuracy v. the test retest scores

%Find common methods
[methods, indsPTR, indsP] = intersect(PTRmethods, Pmethods);

%Loop over methods
H = [ ];
for n = 1:numel(methods)
    
    %If they have the same number of results do the scatter plot
    if numel(PTRresults{indsPTR(n)}) == numel(Presults{indsP(n)})
        
        h = figure; H = [H h]; hold on
        %Plot a line indicating equality
        plot([0 1], [0 1], '--', 'Color', 0.8*[1 1 1])
        %Also plot a curve showing how scores would be related if the inferred
        %maps were just the truth with additive, independent noise
        x = linspace(0,1,250);
        y = sqrt(x);
        plot(x, y, 'r--')
        
        %Now sort the results based on the test-retest scores
        %This is to help match the accuracy scores (there are two for every
        %split-half score)
        PTR = PTRresults{indsPTR(n)}(:); P = Presults{indsP(n)}(:);
        [PTR,i] = sort(PTR); PTR = PTR(1:2:end);
        %Plot the average accuracy for each split-half score
        P = P(i); P = (P(1:2:end)+P(2:2:end))/2;
        plot(PTR, P, 'r.', 'MarkerSize', 10);
        
        %And now plot all the data
        plot(PTRresults{indsPTR(n)}(:), Presults{indsP(n)}(:), ...
            'b.', 'MarkerSize', 10);
        
        %Set axes lims
        xlim([0 1]); ylim([0 1]);
        axis square
        %Labels
        xlabel('Test-retest reliability')
        ylabel('Ground truth accuracy')
        title(methods{n}, 'interpreter', 'none')
        
    end
end

end