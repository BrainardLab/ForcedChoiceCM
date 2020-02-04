function OLRayleighMatchLightSettings(p1, p2, test)
% Compute default lights for OLRayleighMatch experiment 
%% Parameters
% Length of scaling factor adjustment arrays
primaries_length = 201;
test_length = 201;

% Scale primary mixture down from max available by
% this factor.
primaryScaleFactor = 1;
if (primaryScaleFactor > 1)
    error('Do not set primaryScaleFactor greater than 1')
end

% Scaling factor for white light
whiteScale = 0.05;

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;

% Scaling factors for intensity of primary and test lights. The two
% primaries are each scaled from 0 to 1, and the arrays are set up so the
% two scaling factors will always add up to 1. We start with primary2
% scaled up to 1 and primary1 scaled down to 0.
p1Scales = linspace(0,1,primaries_length);
p2Scales = 1-p1Scales;
testScales = linspace(0,1,test_length);
%% Get the OneLight calibration structure
%
% Use it to set some device related parameters.
% numCols is number of mirrors in each OL column
cal = OLGetCalibrationStructure;
[spdLength,settingsLength] = size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors;

%% Initialize arrays for storing precomputed spectra
% These will be cycled through in the adjustments.
% Start-stop arrays hold start positions in the first column and stop
% postions in the second column.
testSpdsNominal = zeros(spdLength,test_length);
testSpdsPredicted = zeros(spdLength,test_length);
testSettings = zeros(settingsLength,test_length);
testStartStops = zeros(test_length,2,numCols);

primarySpdsNominal = zeros(spdLength, primaries_length);
primarySpdsPredicted = zeros(spdLength,primaries_length);
primarySettings = zeros(settingsLength, primaries_length);
primaryStartStops = zeros(primaries_length,2,numCols);

%% Specify and initialize "white" light
%
% This is displayed between the times when subject is comparing
% the primaries and test.
%
% The easiest way to set up something reasonable is to turn
% the OneLight on to half of its max.  We could get fancier
% later and explicitly provide a spectrum.
whitePrimaries = whiteScale * ones(settingsLength, 1);
whiteSpdNominal = OLPrimaryToSpd(cal, whitePrimaries);
whiteSettings = OLPrimaryToSettings(cal, whitePrimaries);
[whiteStarts, whiteStops] = OLSettingsToStartsStops(cal, whiteSettings);

%% Get the "dark" spd, that which comes out when primarys are at 0
darkSpd = OLPrimaryToSpd(cal,zeros(settingsLength, 1));

%% Define desired spectrum for test and the two primary lights
%
% Our code above defines the spectral shape of these.  Here we
% figure out what the most intense version of that shape is,
% that we can get on the OL.
testIncrSpdRel = OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax);
testIncrSpd = OLMaximizeSpd(cal,testIncrSpdRel,'lambda',lambda);
primary1IncrSpdRel = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax);
primary1IncrSpd = OLMaximizeSpd(cal, primary1IncrSpdRel,'lambda',lambda);
primary2IncrSpdRel = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax)/3;
primary2IncrSpd = OLMaximizeSpd(cal, primary2IncrSpdRel,'lambda',lambda);

%% Get set of test lights, varying in intensity
for i = 1:test_length
    testSpdsNominal(:,i) = (testScales(i) * testIncrSpd) + darkSpd;
    [testSettings(:,i),~,testSpdsPredicted(:,i)] = OLSpdToSettings(cal, testSpdsNominal(:,i), 'lambda', lambda);
    [testStartStops(i,1,:),testStartStops(i,2,:)] = OLSettingsToStartsStops(cal, testSettings(:,i));
end

%% Get set of primary lights, varying in relative contribution of two primaries
for i = 1:primaries_length
    primarySpdsNominal(:,i) = (primaryScaleFactor*(p1Scales(i) * primary1IncrSpd) + (p2Scales(i) * primary2IncrSpd)) + darkSpd;
    [primarySettings(:,i),~,primarySpdsPredicted(:,i)] = OLSpdToSettings(cal, primarySpdsNominal(:,i), 'lambda', lambda);
    [primaryStartStops(i,1,:), primaryStartStops(i,2,:)] = OLSettingsToStartsStops(cal, primarySettings(:,i));
end

%% Take a look at spectra (optional)
makePlots = false; 
if makePlots
    figure; clf; hold on; title('Test');
    OLplotSpdCheck(cal.computed.pr650Wls, testSpdsNominal);
    
    figure; clf; title('Primaries');
    OLplotSpdCheck(cal.computed.pr650Wls, primarySpdsNominal);
end
file = sprintf('OLRayleighMatchFineSpectralSettings'); 
% file = sprintf('OLRayleighMatchSpectralSettings_%g_%g_%g.mat', p1, p2, test);
save(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', file));
end 