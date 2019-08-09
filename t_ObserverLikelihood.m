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
wls = SToWls(S);

% Apparatus parameters.  These match the Nagel in wavelengths.
testWavelength = 520;
primaryWavelength1 = 430;
primaryWavelength2 = 545;
primaryWavelength3 = 670;

% Compute indices so that we can set spectra below
testIndex = find(wls == testWavelength);
primaryIndex1 = find(wls == primaryWavelength1);
primaryIndex2 = find(wls == primaryWavelength2);
primaryIndex3 = find(wls == primaryWavelength3);
unitTestSpectrum = zeros(size(wls)); unitTestSpectrum(testIndex) = 1;
unitPrimarySpectrum1 = zeros(size(wls)); unitPrimarySpectrum1(primaryIndex1) = 1;
unitPrimarySpectrum2 = zeros(size(wls)); unitPrimarySpectrum2(primaryIndex2) = 1;
unitPrimarySpectrum3 = zeros(size(wls)); unitPrimarySpectrum3(primaryIndex3) = 1;
primaryBasis = [unitPrimarySpectrum1 unitPrimarySpectrum2 unitPrimarySpectrum3];

% Baseline cone fundamental information
fieldSizeDegrees = 10;
ageInYears = 32;
pupilDiameterMM = 3;

% Individual observer cone fundamentals parameters
observerParams1 = ObserverVecToParams;

% Set color difference model parameters
observerParams1.colorDiffParams.type = 'opponentContrast';
observerParams1.colorDiffParams.LMRatio = 2;
observerParams1.colorDiffParams.lumWeight = 4;
observerParams1.colorDiffParams.rgWeight = 2;
observerParams1.colorDiffParams.byWeight = 0.5;
observerParams1.colorDiffParams.noiseSd = 0.02;

%% Get cone fundamentals, just to make sure they look sensible
[~,T] = ComputeChoiceLikelihood(observerParams1,S,fieldSizeDegrees,ageInYears,pupilDiameterMM);
coneFunadmentalFig = figure; clf; hold on;
plot(SToWls(S),T(1,:)','r','LineWidth',2);
plot(SToWls(S),T(2,:)','g','LineWidth',2);
plot(SToWls(S),T(3,:)','b','LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Fundamental');

%% Find exact metarmer, so that we are working in a reasonable range
M_PrimaryToLMS = T*primaryBasis;
M_LMSToPrimary = inv(M_PrimaryToLMS);
testLMS = T*unitTestSpectrum;
metamerPrimaryWeights = M_LMSToPrimary*testLMS;
metamerSpectrum = primaryBasis*metamerPrimaryWeights;
metamerLMS = T*metamerSpectrum;
if (max(abs(metamerLMS-testLMS)./testLMS) > 1e-6)
    error('Failed to compute good metamer');
end

%% Let's compute some likelihoods
%
% We keep the first comparison a match and look at the probability it
% is chosen as we march the other other by incrementing one of the
% primary values.  Expect probs1 to start at 0.5 and increase to 1.
reference = unitTestSpectrum;
primaryDeltaFactor = 1/200;
deltaPrimary = primaryDeltaFactor*metamerPrimaryWeights(1);
comparison1 = primaryBasis*metamerPrimaryWeights;
nSteps = 100;
for ii = 0:nSteps-1
    comparison2PrimaryWeights = metamerPrimaryWeights + [ii*deltaPrimary 0 0]';
    comparison2 = primaryBasis*comparison2PrimaryWeights;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams1,S,fieldSizeDegrees,ageInYears,pupilDiameterMM,...
        reference,comparison1,comparison2);
end
likelihoodFigure2 = figure; clf; hold on
plot(1:nSteps,probs1,'ro','MarkerSize',12','MarkerFaceColor','r');
xlabel('Step Size');
ylabel('Probability 1 Chosen');

% We keep push the first comparison away from the test, and march the
% second from the test through the first and beyond. Expect probs1 to
% start at 0 and increase to 1, with 50% point halfway.
reference = unitTestSpectrum;
primaryDeltaFactor = 1/200;
deltaPrimary = primaryDeltaFactor*metamerPrimaryWeights(1);
comparison1 = primaryBasis*(metamerPrimaryWeights + [nSteps/2*deltaPrimary 0 0]');
nSteps = 100;
for ii = 0:nSteps-1
    comparison2PrimaryWeights = metamerPrimaryWeights + [ii*deltaPrimary 0 0]';
    comparison2 = primaryBasis*comparison2PrimaryWeights;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams1,S,fieldSizeDegrees,ageInYears,pupilDiameterMM,...
        reference,comparison1,comparison2);
end
likelihoodFigure2 = figure; clf; hold on
plot(1:nSteps,probs1,'ro','MarkerSize',12','MarkerFaceColor','r');
xlabel('Step Size');
ylabel('Probability 1 Chosen');

% Now run this with a different observer. There should be some effect,
% and there is.  It is a bit harder for me to intuit what the effect 
% should be, so we'll just be happy that there is some effect.
observerParams2 = observerParams1;
observerParams2.indDiffParams.lambdaMaxShift(1) = 5;
observerParams2.indDiffParams.lambdaMaxShift(2) = -5;
observerParams2.indDiffParams.lambdaMaxShift(3) = -10;
reference = unitTestSpectrum;
primaryDeltaFactor = 1/200;
deltaPrimary = primaryDeltaFactor*metamerPrimaryWeights(1);
comparison1 = primaryBasis*(metamerPrimaryWeights + [nSteps/2*deltaPrimary 0 0]');
nSteps = 100;
for ii = 0:nSteps-1
    comparison2PrimaryWeights = metamerPrimaryWeights + [ii*deltaPrimary 0 0]';
    comparison2 = primaryBasis*comparison2PrimaryWeights;
    probs1(ii+1) = ComputeChoiceLikelihood(observerParams2,S,fieldSizeDegrees,ageInYears,pupilDiameterMM,...
        reference,comparison1,comparison2);
end
figure(likelihoodFigure2);
plot(1:nSteps,probs1,'bo','MarkerSize',12','MarkerFaceColor','b');







