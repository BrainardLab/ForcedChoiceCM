% Illustrate and test likelihood computation routines
%
% Description:
%     Illustrates how we use the routines here to compute
%     likelihood in a color selection version of
%     the fundamental color matching experiment

% History:
%   08/08/19  dhb  Wrote it.

%% Clear
clear; close all;

%% Set up parameters

% Apparatus parameters
apparatusParams = DefaultApparatusParams('monochromatic');

% Individual observer cone fundamentals parameters
observerParams1.coneParams = DefaultConeParams('cie_asano');

% Set color difference model parameters
observerParams1.colorDiffParams = DefaultColorDiffParams('opponentContrast');

%% Get cone funadmentals
T = ComputeObserverFundamentals(observerParams1.coneParams,apparatusParams.S);

%% Find exact metamer, so that we are working in a reasonable range
adaptation = apparatusParams.unitTestSpectrum;
reference = apparatusParams.unitTestSpectrum;
[metamerPrimary,metamerSpectrum,metamerLMS] = FindMetamer(apparatusParams,T,reference);

%% Let's compute some likelihoods
%
% We keep the first comparison a match and look at the probability it
% is chosen as we march the other other by incrementing one of the
% primary values.  Expect probs1 to start at 0.5 and increase to 1.
primaryDeltaFactor = 1/200;
deltaPrimary = primaryDeltaFactor*metamerPrimary(1);
comparison1 = apparatusParams.primaryBasis*metamerPrimary;
nSteps = 100;
for ii = 0:nSteps-1
    comparison2PrimaryWeights = metamerPrimary + [ii*deltaPrimary 0 0]';
    comparison2 = apparatusParams.primaryBasis*comparison2PrimaryWeights;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams1,apparatusParams.S,...
        adaptation,reference,comparison1,comparison2);
end
likelihoodFigure2 = figure; clf; hold on
plot(1:nSteps,probs1,'ro','MarkerSize',12','MarkerFaceColor','r');
xlabel('Step Size');
ylabel('Probability 1 Chosen');

% We keep push the first comparison away from the reference, and march the
% second from the test through the first and beyond. Expect probs1 to
% start at 0 and increase to 1, with 50% point halfway.
primaryDeltaFactor = 1/200;
deltaPrimary = primaryDeltaFactor*metamerPrimary(1);
comparison1 = apparatusParams.primaryBasis*(metamerPrimary + [nSteps/2*deltaPrimary 0 0]');
nSteps = 100;
for ii = 0:nSteps-1
    comparison2PrimaryWeights = metamerPrimary + [ii*deltaPrimary 0 0]';
    comparison2 = apparatusParams.primaryBasis*comparison2PrimaryWeights;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams1,apparatusParams.S,...
        adaptation,reference,comparison1,comparison2);
end
likelihoodFigure2 = figure; clf; hold on
plot(1:nSteps,probs1,'ro','MarkerSize',12','MarkerFaceColor','r');
xlabel('Step Size');
ylabel('Probability 1 Chosen');

% Now run this with a different observer. There should be some effect,
% and there is.  It is a bit harder for me to intuit what the effect 
% should be, so we'll just be happy that there is some effect.
observerParams2 = observerParams1;
observerParams2.coneParams.indDiffParams.lambdaMaxShift(1) = 5;
observerParams2.coneParams.indDiffParams.lambdaMaxShift(2) = -5;
observerParams2.coneParams.indDiffParams.lambdaMaxShift(3) = -10;
primaryDeltaFactor = 1/200;
deltaPrimary = primaryDeltaFactor*metamerPrimary(1);
comparison1 = apparatusParams.primaryBasis*(metamerPrimary + [nSteps/2*deltaPrimary 0 0]');
nSteps = 100;
for ii = 0:nSteps-1
    comparison2PrimaryWeights = metamerPrimary + [ii*deltaPrimary 0 0]';
    comparison2 = apparatusParams.primaryBasis*comparison2PrimaryWeights;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams2,apparatusParams.S,...
        adaptation,reference,comparison1,comparison2);
end
figure(likelihoodFigure2);
plot(1:nSteps,probs1,'bo','MarkerSize',12','MarkerFaceColor','b');







