% t_AnalyticRayleighMatch
%
% Tutorial/test program to illustrate analytic way of computing Rayleigh
% match and compare with other ways of doing so.
%
% 08/07/2020  dhb  Wrote it.

%% Clear
clear; close all;

%% Define primary and test peak wavelengths
p1Wl = 670;
p2Wl = 560;
testWl = 600;
p1ScaleFactor = 1;
p2ScaleFactor = 0.02;
testScaleFactor = 0.5;
adjustmentLength = 201;

%% Get OL spectra
startsStopsFilename = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
    adjustmentLength, p1Wl, p2Wl, testWl, p1ScaleFactor, p2ScaleFactor, testScaleFactor);

% Precompute if needed
if (~exist(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'precomputedStartStops', startsStopsFilename)))
    fprintf('Precomputing starts/stops\n');
    OLRayleighMatchLightSettings(p1Wl, p2Wl, testWl, ...
        'p1ScaleFactor',p1ScaleFactor,'p2ScaleFactor',p2ScaleFactor,'testScaleFactor',testScaleFactor, ...
        'adjustmentLength',adjustmentLength);
end

% Load the various OL spectra and related info.
fprintf('Loading precomputed starts/stops\n');
olData = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'precomputedStartStops', startsStopsFilename));
S = WlsToS(olData.wls);

%% Some plots of spectra

% What we're plotting
%   olData.primarySpdsPredicted - Mixtures of primaries available in the experiment
%   olData.primary1IncrSpd      - Nominal primary 1 at full power, dark light subtracted
%   olData.primary2IncrSpd      - Nominal primary 2 at full power, dark light subtracted
figure; clf; hold on
plot(olData.wls,olData.primarySpdsPredicted-olData.darkSpd);
plot(olData.wls,p1ScaleFactor*olData.primary1IncrSpdPredicted,'r','LineWidth',2);
plot(olData.wls,p2ScaleFactor*olData.primary2IncrSpdPredicted,'g','LineWidth',2);
plot(olData.wls,p1ScaleFactor*olData.primary1IncrSpd,'k:','LineWidth',1);
plot(olData.wls,p2ScaleFactor*olData.primary2IncrSpd,'k:','LineWidth',1);

% What we're plotting
%    olData.testSpdsPredicted  - Test lights available in the experiment
%    olData.testIncrSpd        - Target test light at full power, dark light subtracted
figure; clf; hold on
plot(olData.wls,olData.testSpdsPredicted-olData.darkSpd);
plot(olData.wls,testScaleFactor*olData.testIncrSpdPredicted,'y','LineWidth',2);
plot(olData.wls,testScaleFactor*olData.testIncrSpd,'k:','LineWidth',2);

%% Get an observer
observer = genRayleighObserver('S',S);

%% Compute analytic predicted Rayleigh match
[testAdjustedSpd,primaryMixtureSpd,testIntensity,lambda] = ...
    computePredictedRayleighMatch(p1Spd,p2Spd,testSpd,ObserverParamsToVec(observer.coneParams.indDiffParams),opponentParams, ...
    'fieldSize',observer.coneParams.fieldSize,'age',observer.coneParams.age,'S',S);
