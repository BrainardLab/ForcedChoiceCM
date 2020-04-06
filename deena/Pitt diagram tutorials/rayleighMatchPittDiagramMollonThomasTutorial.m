% Illustrates how Rayleigh matches vary with changes in cone lambda max and
% optical density
%
% Description:
%   Simulate Rayleigh match performance for variations in L and M cone
%   peak sensitivity and optical density, as described by Thomas and Mollon
%   (2003). Plots Rayleigh match performance in the form of a Pitt diagram.
%
%   The simulated anomaloscope allows adjustment of a monochromatic test
%   and the ratio of two monochromatic primaries in the match. The routine
%   computes the cone responses to the test and match and from these a
%   color difference.  Matches are predicted for test intensity and mixing
%   ratio parameters where the color difference is below a criterion. The
%   locus of matches is plotted in a Pitt diagram, where the x-axis is the
%   mixing ratio and the y-axis is the test intensity.
%
%   The routine produces two plots that come close to the Thomas and
%   Mollon results. The first plots Rayleigh match performance for
%   hypothetical dichromats with varying L/M cone lambda maxes, while the
%   second plots Rayleigh matches for dichromats with varying optical
%   densities. The routine also makes two plots that compare fit lines of 
%   the simulated results for lambda max and optical density with fit lines 
%   extracted from theThomas and Mollon plots (using WebPlotDigitizer), as
%   well as two plots that compare the slopes of these fit lines for each 
%   condition.

% History
%   07/03/19  dhb       Wrote it.
%   09/03/19  dhb, dce  Modified to use Asano et al. individual difference
%                       parameters, but in the end does the same thing.
%                       However, and enterprising person can now examine
%                       the effect of changing photopigment density.
%   09/12/19  dce       Modified to include specific plots of Thomas and
%                       Mollon results
%   10/2/19   dce       Added comments and comparative routines
%   10/10/19  dce       Cleaned up. Made ComputeRayleighConfusion a
%                       separate routine
%% Clear
clear; close all;

%% In our first case, we will plots effects of variation in cone lambda max
%  for a hypothetical dichromat. We begin by defining our parameters

% Cone lambda maxes, taken from Fig. 1 of the Thomas and Mollon paper.
lambdaMaxes1 = [  ...
    [531 531 420.7]' ...
    [541 541 420.7]' ...
    [551 551 420.7]' ...
    [561 561 420.7]'];

% We specify the cones as a shift relative to a nomogram generated lambda
%  max. These base values are set equal to the test values
baseLambdaMaxes1 = lambdaMaxes1;

% You can also allow the specified photopigment density to vary. Enter
% these as percent changes relative to nominal values (positive or
% negative). Here, all photopigment densities are set to zero.
dphotopigments1 = zeros(3,4);

% Define plotting parameters and produce Figure 1 
colors1 = [ 'r' 'g' 'b' 'y'];
legend1 = {'531' '541' '551' '561'};
title1 = 'Lambda Max Variation';
plotMatches(lambdaMaxes1, baseLambdaMaxes1, dphotopigments1, colors1,...
    legend1, title1, true)

%% The second case plots effects of optical density variation for a
%% hypothetical dichromat

% All lambda max and base lambda max values are the same
lambdaMaxes2 = [[561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]'];
baseLambdaMaxes2 = lambdaMaxes2;

% Optical density shifts are taken from Fig. 1 of the Thomas and Mollon
% paper: 0.05, 0.2, 0.35, 0.5, 0.65. The values here are calculated as
% percent changes from the average optical density of 0.3.
dphotopigments2 = [  ...
    [-83.33 -83.33 0]' ...
    [-33.33 -33.33 0]' ...
    [16.67 16.67 0]' ...
    [66.67 66.67 0]'...
    [116.67 116.67 0]'];

% Define plotting parameters and plot Figure 2
colors2 = [ 'r' 'g' 'b' 'y' 'm' ];
legend2 = { '-83.33%' '-33.33%' '16.67%' '66.67%' '116.67%' };
title2 = 'Optical Density Variation';
plotMatches(lambdaMaxes2, baseLambdaMaxes2, dphotopigments2, colors2,...
    legend2, title2, false);

% Finish off the Rayleigh match comparison plots (Figures 3-4)
for i = 3:4
    figure(i);
    h = findobj(figure(i),'Type','line');
    last = length(h);
    legend([h(last) h(1)], {'Mollon Thomas Matches','Simulated Matches'}); 
end 

% Finish off the slope comparison plots (Figures 5-6)
figure(5);
hold on;
plot(-0.4:0.01:0.2, -0.4:0.01:0.2, 'r-'); % comparison 1:1 line  
title('Slope comparison - lambda max variation');
xlabel('Simulated Rayleigh match slope');
ylabel('Mollon-Thomas Rayleigh match slope');
legend(legend1)

figure(6);
hold on; 
plot(-0.2:0.01:0.2, -0.2:0.01:0.2, 'r-'); % comparison 1:1 line  
title('Slope comparison - optical density variation');
xlabel('Simulated Rayleigh match slope');
ylabel('Mollon-Thomas Rayleigh match slope'); 
legend(legend2)

function plotMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments,...
    colors, theLegend, theTitle, lambdaMax)
% Create Pitt diagrams based on passed cone parameters
% Syntax:
%   plotMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments, colors,...
%                theLegend, the Title, od)
% Description:
%    Takes in cone and plot parameters for a series of observers, then uses
%    the anomaloscope simulation code to produce a Pitt diagram overlaying
%    the regions of confusion for each observer.
%
% Inputs
%    lambdaMaxes      - 3xn array of peak cone spectral sensitivities (nm)
%                       for each of the n simulated observers you will be
%                       plotting. The first row contains l cone values, the
%                       second m cone values, and the third s cone values.
%    baseLambdaMaxes  - comparison array for lambdaMaxes. Usually set
%                       equal to lambdaMaxes
%    dphotopigments   - 3xn array of optical density variations, expressed
%                       as percent changes from an average value. The first
%                       row contains l cone values, the second contains m
%                       cone values, and the third contains s cone values
%    colors           - 1xn character vector containing keys for plot colors
%    theLegend        - 1x1 cell array containing n character vectors for
%                       the plot legend entries.
%    theTitle         - character vector for lot title
%    lambdaMax        - logical indicating whether we are plotting lambda
%                       max variation (true) or optical density variation
%                       (false).
%
% Outputs:
%    Four plots, as described at the top of the file

% Range parameters for simulated anomaloscope
% Mixing ratio of zero is all green primary, 1 is all red
% (assuming shorter primary is first one specified in routine below.)
mixingRatioRange = 0:0.001:1;
testIntensityRange = 0:0.001:0.35;
testWavelength = 590; 
thresholdVal = 0.12; % Threshold difference below which it is a match

% Find the individual different parameters for each observer
for ii = 1:size(lambdaMaxes,2)
    indDiffParams(ii).dlens= 0;
    indDiffParams(ii).dmac = 0;
    indDiffParams(ii).dphotopigment = dphotopigments(:,ii)';
    indDiffParams(ii).lambdaMaxShift = lambdaMaxes(:,ii)' - baseLambdaMaxes(:,ii)';
    indDiffParams(ii).shiftType = 'linear';
end

% Initialize plots comparing simulation to Mollon-Thomas results, and plot
% the Mollon-Thomas data in preparation
if lambdaMax % Lambda max variation
    load mollon_thomas_lambdaMaxFitParams;
    fitParams = lambdaMaxFitParams;
    figure(3);
    title('Rayleigh match variation with lambda max');
else % Optical density variation
    load mollon_thomas_odFitParams;
    fitParams = odFitParams;
    figure(4);
    title('Rayleigh match variation with optical density density');
end
% Plot Mollon-Thomas data
hold on;
[row, ~] = size(fitParams);
for i = 1:row
    plot(mixingRatioRange, fitParams(i,1) * mixingRatioRange...
        + fitParams(i,2), 'b');
end
xlabel('Mixing Ratio (0 -> green; 1 -> red)');
ylabel('Test Intensity');
hold off;

% Plot simulated Rayleigh matches for each observer by using the
% ComputeConfusions routine and the passed parameters
theFigure = figure; clf; hold on
for kk = 1:size(lambdaMaxes,2)
    [testIntensity{kk},mixingRatio{kk},matchDiff{kk}] =...
        ComputeRayleighConfusions(baseLambdaMaxes(:,kk),indDiffParams(kk),...
        testWavelength,testIntensityRange,mixingRatioRange);
    figure(theFigure);
    index = find(matchDiff{kk} < thresholdVal);
    plot(mixingRatio{kk}(index),testIntensity{kk}(index),...
        [colors(kk) 'o'],'MarkerFaceColor',colors(kk));
    
    % Add these simulated Rayleigh matches to the plot comparing simulated
    % and Mollon-Thomas results (+0.062 for od, + 0.08 for lambda max)
    if lambdaMax; figure(3); else; figure(4); end
    hold on; 
    tutorialFit = polyfit(mixingRatio{kk}(index),...
        testIntensity{kk}(index), 1);
    plot(mixingRatioRange, mixingRatioRange * tutorialFit(1)...
        + tutorialFit(2), 'r');
    hold off;
    
    % Compare the slopes of simulated and Mollon-Thomas matches for each
    % observer
    if lambdaMax; figure(5); else; figure(6); end 
    hold on;
    plot(tutorialFit(1), fitParams(kk, 1), 'Color', colors(kk), 'LineStyle', 'none', 'Marker', '*');
    hold off; 
end

% Finish off the Pitt diagram plots (Figures 1-2)
figure(theFigure);
xlim([min(mixingRatioRange) max(mixingRatioRange)]);
ylim([min(testIntensityRange) max(testIntensityRange)]);
xlabel(' Mixing Ratio (0 -> green; 1 -> red)');
ylabel('Test Intensity');
axis('square');
legend(theLegend);
title(theTitle);
end
