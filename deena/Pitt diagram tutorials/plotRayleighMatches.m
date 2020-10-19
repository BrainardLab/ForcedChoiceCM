function plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments,...
    colors, theLegend, theTitle, varargin)

% Create Pitt diagrams based on passed cone parameters
%
% Syntax:
%   plotMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments, colors,...
%                theLegend, the Title)
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
% Outputs:
%    Pitt Diagram that plots passed cone parameters
%
% Optional key-value pairs:
%    'test'         - integer wavelength of the test light in nm. Default
%                     is 590.
%    'p1'           - integer wavelength of the first primary light in nm.
%                     Default is 545.
%    'p2'           - integer wavelength of the second primary light in nm.
%                     Default is 679.

% History
%    xx/xx/19   dce  - Wrote program
%    4/5/20     dce  - Edited for style

% Parse input
p = inputParser;
p.addParameter('p1', 560, @(x) (isnumeric(x)));
p.addParameter('p2', 670, @(x) (isnumeric(x)));
p.addParameter('test', 590, @(x) (isnumeric(x)));
p.parse(varargin{:});

% Range parameters for simulated anomaloscope
% Mixing ratio of zero is all short primary primary, 1 is all long
% (assuming shorter primary is first one specified in routine below.)
mixingRatioRange = 0:0.001:1;
testIntensityRange = 0:0.001:0.35;
thresholdVal = 0.12; % Threshold difference below which it is a match

% Find the individual different parameters for each observer
for ii = 1:size(lambdaMaxes,2)
    indDiffParams(ii).dlens= 0;
    indDiffParams(ii).dmac = 0;
    indDiffParams(ii).dphotopigment = dphotopigments(:,ii)';
    indDiffParams(ii).lambdaMaxShift = lambdaMaxes(:,ii)' - baseLambdaMaxes(:,ii)';
    indDiffParams(ii).shiftType = 'linear';
end

% Plot simulated Rayleigh matches for each observer by using the
% ComputeConfusions routine and the passed parameters
theFigure = figure; clf; hold on
for kk = 1:size(lambdaMaxes,2)
    [testIntensity{kk},mixingRatio{kk},matchDiff{kk}] = ComputeRayleighConfusions(baseLambdaMaxes(:,kk),indDiffParams(kk),p.Results.p1, p.Results.p2, p.Results.test, testIntensityRange,mixingRatioRange);
    figure(theFigure);
    index = find(matchDiff{kk} < thresholdVal);
    plot(mixingRatio{kk}(index),testIntensity{kk}(index),...
        [colors(kk) 'o'],'MarkerFaceColor',colors(kk));
    
    % Finish off the Pitt diagram plots (Figures 1-2)
    figure(theFigure);
    xlim([min(mixingRatioRange) max(mixingRatioRange)]);
    ylim([min(testIntensityRange) max(testIntensityRange)]);
    xlabel(' Primary Ratio (Proportion Red)');
    ylabel('Test Intensity');
    axis('square');
    legend(theLegend);
    title(theTitle);
end
end