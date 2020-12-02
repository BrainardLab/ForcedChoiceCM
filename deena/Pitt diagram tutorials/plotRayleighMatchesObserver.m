function plotRayleighMatchesObserver(observer,p1Spd,p2Spd,testSpd,...
    noiseScalar,colorName,theTitle,varargin)
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
%                       (green)
%    p2Spd            - Column vector representation of the second primary
%                       (red)
%    testSpd          - Column vector representation of the test spd 
%    colorName        - Character vector with color name
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
%    11/18/20   dce  - Edited plotting

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
hold on;

[testIntensity,mixingRatio,matchDiff] = ...
    ComputeRayleighConfusionsObs(observer,p1Spd,p2Spd,testSpd,...
    testIntensityRange,mixingRatioRange);

index = find(matchDiff < thresholdVal);
plot(mixingRatio(index),testIntensity(index),'.','MarkerSize',20,...
    'MarkerFaceColor',colorName,'MarkerEdgeColor',colorName);

% Finish off the Pitt diagram plots (Figures 1-2)
xlim([0 1]);
ylim([0 1]);
xlabel(' Primary Ratio (Proportion Red)');
ylabel('Test Intensity');
axis('square');
title(theTitle,'interpreter','Latex');
end