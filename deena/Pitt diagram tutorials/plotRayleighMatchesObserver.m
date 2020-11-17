function plotRayleighMatchesObserver(observer,p1Spd,p2Spd,testSpd,...
    noiseScalar,colors,theLegend,theTitle,varargin)
% Create Pitt diagrams using parameters from the observer model
%
% Syntax:
%   plotMatches(observer,colors,theLegend, the Title)
% Description:
%    Produces a Pitt diagram using the observer model. Takes in an observer
%    struct, Rayleigh matching spectra, and plotting parameters (colors, 
%    legend, and title). Uses the computeRayleighConfusionsObserver routine
%    to compute a Pitt diagram for regions of observer confusion.
%
% Inputs
%    observer         - Struct with observer parameters. Must contain the
%                       following fields: colorDiffParams, T_cones. See
%                       genRayleighObserver for details.
%    p1Spd            - Column vector representation of the first primary
%    p2Spd            - Column vector representation of the second primary
%    testSpd          - Column vector representation of the test spd 
%    colors           - 1xn character vector containing keys for plot colors
%    theLegend        - 1x1 cell array containing n character vectors for
%                       the plot legend entries.
%    theTitle         - Character vector for plot title
%
% Outputs:
%    None (produces a figure).
%
% Optional key-value pairs:
%    'figHandle'     - Optional figure handle. If nonempty, plot is made in
%                      the figure handle instead of in a new figure. 
%                      Default is [].

% History
%    11/16/20   dce  - Adapeed from plotRayleighMatches

% Parse input
p = inputParser;
p.addParameter('figHandle', []);
p.parse(varargin{:});

% Range parameters for simulated anomaloscope
% Mixing ratio of zero is all short primary primary, 1 is all long
% (assuming shorter primary is first one specified in routine below.)
mixingRatioRange = 0:0.01:1;
testIntensityRange = 0:0.01:1;
thresholdVal = observer.colorDiffParams.noiseSd * noiseScalar; 

% Plot simulated Rayleigh matches for each observer by using the
% ComputeConfusions routine and the passed parameters
if ~isempty(p.Results.figHandle)
    theFigure = p.Results.figHandle;
else 
    theFigure = figure();
end 
clf; hold on

[testIntensity,mixingRatio,matchDiff] = ...
    ComputeRayleighConfusionsObs(observer,p1Spd,p2Spd,testSpd,...
    testIntensityRange,mixingRatioRange);
% Reshape arrays 
testIntensity = reshape(testIntensity,[1 length(mixingRatioRange)...
    *length(testIntensityRange)]);
mixingRatio = reshape(mixingRatio,[1 length(mixingRatioRange)...
    *length(testIntensityRange)]);
matchDiff = reshape(matchDiff,[1 length(mixingRatioRange)...
    *length(testIntensityRange)]);

index = find(matchDiff < thresholdVal);
plot(mixingRatio(index),testIntensity(index),...
    [colors 'o'],'MarkerFaceColor',colors);

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