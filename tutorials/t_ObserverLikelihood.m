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
%
% Wavelength sampling
S = [400 1 301];

% Apparatus parameters
stimulusParams.matchApparatusParams = DefaultMatchApparatusParams('monochromatic',S);
stimulusParams.testParams = DefaultTestParams('monochromatic',S);

% Individual observer cone fundamentals parameters
observerParams1.coneParams = DefaultConeParams('cie_asano');

% Set color difference model parameters
observerParams1.colorDiffParams = DefaultColorDiffParams('opponentContrast');

%% Get cone funadmentals
T = ComputeObserverFundamentals(observerParams1.coneParams,S);

%% Find exact metamer, so that we are working in a reasonable range
adaptationSpd = stimulusParams.testParams.testIntensity*stimulusParams.testParams.unitTestSpectrum;
referenceSpd = stimulusParams.testParams.testIntensity*stimulusParams.testParams.unitTestSpectrum;
[metamerPrimary,metamerSpd] = FindMetamer(stimulusParams.matchApparatusParams,T,T*referenceSpd);
metamerLMS = T*metamerSpd;
comparison1Spd = stimulusParams.matchApparatusParams.primaryBasis*metamerPrimary;

%% Find comparisons in interesting directions and sizes
%
% This is one noise sd worth of step in the r/g direction.
opponentDirContrast = [0,observerParams1.colorDiffParams.noiseSd 0]';
opponentDirLMS = OpponentContrastToLMS(GetOpponentContrastMatrix(observerParams1.colorDiffParams),T*metamerSpd,opponentDirContrast);
[opponentDirPrimary,opponentDirSpd] = FindMetamer(stimulusParams.matchApparatusParams,T,opponentDirLMS);
opponentDirLMSCheck = T*opponentDirSpd;
if (max(abs(opponentDirLMS-opponentDirLMSCheck)) > 1e-8)
    error('Do not recover metameric cone excitations properly');
end
[opponentDirDiff,opponentDirContrastCheck] = ComputeMatchDiff(observerParams1.colorDiffParams, ...
    GetOpponentContrastMatrix(observerParams1.colorDiffParams),T*adaptationSpd,T*referenceSpd,opponentDirLMS);
if (max(abs(opponentDirContrast-opponentDirContrastCheck)) > 1e-8)
    error('Do not recover opponent direction contrast the way we should');
end
opponentDirSpdCheck = stimulusParams.matchApparatusParams.primaryBasis*opponentDirPrimary;
if (max(abs(opponentDirSpd - opponentDirSpdCheck)) > 1e-8)
    error('Do not recover opponent direction spectrum properly');
end
opponentDirLMSCheck2 = T*opponentDirSpdCheck;
if (max(abs(opponentDirLMS-opponentDirLMSCheck2)) > 1e-8)
    error('Do not recover opponent direction cone excitations properly');
end
opponentDirDeltaPrimary = opponentDirPrimary-metamerPrimary;
opponentDirDeltaSpd = opponentDirSpd-metamerSpd;
probCheck = ComputeChoiceLikelihood(observerParams1, GetOpponentContrastMatrix(observerParams1.colorDiffParams), ...
        T*adaptationSpd,T*referenceSpd,T*comparison1Spd,T*(comparison1Spd + opponentDirDeltaSpd));
if (abs(probCheck - normcdf(1,0,1)) > 1e-8)
    error('Did not properly construct a one noise sd opponent step.');
end

%% Let's compute some likelihoods
%
% We keep the first comparison a match and look at the probability it
% is chosen as we march the other other by incrementing one of the
% primary values.  Expect probs1 to start at 0.5 and increase to 1.
nSteps = 100;
primaryDeltaFactor = 2/nSteps;   
for ii = 0:nSteps-1
    comparison2Primary = metamerPrimary + ii*primaryDeltaFactor*opponentDirDeltaPrimary;
    comparison2Spd = stimulusParams.matchApparatusParams.primaryBasis*comparison2Primary;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams1, GetOpponentContrastMatrix(observerParams1.colorDiffParams), ...
        T*adaptationSpd,T*referenceSpd,T*comparison1Spd,T*comparison2Spd);
end
likelihoodFigure1 = figure; clf; hold on
plot(1:nSteps,probs1,'ro','MarkerSize',12','MarkerFaceColor','r');
xlabel('Step Size');
ylabel('Probability 1 Chosen');

% We keep push the first comparison away from the reference, and march the
% second from the test through the first and beyond. Expect probs1 to
% start at 0 and increase to 1, with 50% point halfway.
comparison1Spd = stimulusParams.matchApparatusParams.primaryBasis*(metamerPrimary + (nSteps/2)*primaryDeltaFactor*opponentDirDeltaPrimary);
for ii = 0:nSteps-1
    comparison2Primary = metamerPrimary + ii*primaryDeltaFactor*opponentDirDeltaPrimary;
    comparison2Spd = stimulusParams.matchApparatusParams.primaryBasis*comparison2Primary;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams1, GetOpponentContrastMatrix(observerParams1.colorDiffParams), ...
        T*adaptationSpd,T*referenceSpd,T*comparison1Spd,T*comparison2Spd);
end
likelihoodFigure2 = figure; clf; hold on
plot(0:nSteps-1,probs1,'ro','MarkerSize',12','MarkerFaceColor','r');
xlabel('Step Size');
ylabel('Probability 1 Chosen');

% Now run this with a different observer. There should be some effect,
% and there is.  It is a bit harder for me to intuit what the effect 
% should be, so we'll just be happy that there is some effect.
observerParams2 = observerParams1;
observerParams2.coneParams.indDiffParams.lambdaMaxShift(1) = 3;
observerParams2.coneParams.indDiffParams.lambdaMaxShift(2) = 0;
observerParams2.coneParams.indDiffParams.lambdaMaxShift(3) = 0;
for ii = 0:nSteps-1
    comparison2Primary = metamerPrimary + ii*primaryDeltaFactor*opponentDirDeltaPrimary;
    comparison2Spd = stimulusParams.matchApparatusParams.primaryBasis*comparison2Primary;
    T2 = ComputeObserverFundamentals(observerParams2.coneParams,S);
    probs2(ii+1) = ComputeChoiceLikelihood(observerParams2, GetOpponentContrastMatrix(observerParams2.colorDiffParams), ...
        T2*adaptationSpd,T2*referenceSpd,T2*comparison1Spd,T2*comparison2Spd);
end
figure(likelihoodFigure2);
plot(0:nSteps-1,probs2,'bo','MarkerSize',12','MarkerFaceColor','b');







