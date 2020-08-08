function OLRayleighMatchLightSettings(p1, p2, test, varargin)
% Compute spectra to display for OLRayleighMatch experiment
% Syntax:
%   OLRayleighMatch
%
% Description:
%    Given two primary wavelengths and a test wavelength, computes a range
%    of precomputed spectra to be displayed with various mixtures of the
%    primary lights and various intensities of the test lights. These
%    results are saved as OneLight start/stops arrays which can be used to
%    run a Rayleigh matching experiment.
%
%    Scaling factors for the two primaries and the test light can be
%    entered as optional key-value pairs
%
% Inputs:
%    'p1'        - integer wavelength of the first primary light in nm.
%    'p2'        - integer wavelength of the second primary light in nm.
%    'test'      - integer wavelength of the test light in nm.
%
% Outputs:
%    Saves a .mat file which includes computed spds, start/stops data,
%    white light, and other experimental parameters.
%
% Optional key-value pairs
%    'p1ScaleFactor'    - Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2ScaleFactor'    - Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 1.
%    'testScaleFactor'  - Numerical scale factor for the test light,
%                         between 0 and 1. Default is 1.
%    'whiteScaleFactor' - Numerical scale factor for the white light,
%                         between 0 and 1. Default is 0.05.
%    'adjustmentLength' - length of scaling factor adjustment arrays.
%                         Default is 201.

% History:
%   02/4/20    dce       Separated from OLRayleighMatch
%   06/09/20   dce       Added scaling factors as key-value pairs
%   07/14/20   dce       Removed division by 3 from p2 calculation, removed
%                        figure saving
%   08/08/20   dhb       Add predicted incremental primaries and test

%% Parse input
p = inputParser;
p.addParameter('p1ScaleFactor', 1, @(x) (isnumeric(x)));
p.addParameter('p2ScaleFactor', 1, @(x) (isnumeric(x)));
p.addParameter('testScaleFactor', 1, @(x) (isnumeric(x)));
p.addParameter('whiteScaleFactor', 0.05, @(x) (isnumeric(x)));
p.addParameter('adjustmentLength', 201, @(x) (isnumeric(x)));
p.parse(varargin{:});

%% Set up parameters
% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;

adjustmentLength = p.Results.adjustmentLength;
p1ScaleFactor = p.Results.p1ScaleFactor;
p2ScaleFactor = p.Results.p2ScaleFactor;
testScaleFactor = p.Results.testScaleFactor;
whiteScaleFactor = p.Results.whiteScaleFactor;

% Check that scale factor inputs make sense
if (p1ScaleFactor > 1 || p2ScaleFactor > 1 || testScaleFactor > 1 ||...
        whiteScaleFactor > 1)
    error('Scaling factors cannot be greater than 1')
end
if (p1ScaleFactor <= 0 || p2ScaleFactor <= 0 || testScaleFactor <= 0 ||...
        whiteScaleFactor <= 0)
    error('Scaling factors must be greater than 0')
end

% Scaling factors for intensity of primary and test lights. The two
% primaries are each scaled from 0 to 1, and the arrays are set up so the
% two scaling factors will always add up to 1. We start with primary2
% scaled up to 1 and primary1 scaled down to 0.
p1Scales = linspace(0,1,adjustmentLength);
p2Scales = 1-p1Scales;
testScales = linspace(0,1,adjustmentLength);

%% Get the OneLight calibration structure
%
% Use it to set some device related parameters.
% numCols is number of mirrors in each OL column
cal = OLGetCalibrationStructure('CalibrationType',getpref('ForcedChoiceCM','currentCal'));
[spdLength,settingsLength] = size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors;
wls = cal.computed.pr650Wls;

%% Initialize arrays for storing precomputed spectra
% These will be cycled through in the adjustments.
% Start-stop arrays hold start positions in the first column and stop
% postions in the second column.
testSpdsNominal = zeros(spdLength,adjustmentLength);
testSpdsPredicted = zeros(spdLength,adjustmentLength);
testSettings = zeros(settingsLength,adjustmentLength);
testStartStops = zeros(adjustmentLength,2,numCols);

primarySpdsNominal = zeros(spdLength, adjustmentLength);
primarySpdsPredicted = zeros(spdLength,adjustmentLength);
primarySettings = zeros(settingsLength, adjustmentLength);
primaryStartStops = zeros(adjustmentLength,2,numCols);

%% Specify and initialize "white" light
%
% This is displayed between the times when subject is comparing
% the primaries and test.
%
% The easiest way to set up something reasonable is to turn
% the OneLight on to half of its max.  We could get fancier
% later and explicitly provide a spectrum.
whitePrimaries = whiteScaleFactor * ones(settingsLength, 1);
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
[~,~,testIncrSpdPredicted] = OLSpdToSettings(cal,testIncrSpd+darkSpd,'lambda', lambda);
testIncrSpdPredicted = testIncrSpdPredicted-darkSpd;

primary1IncrSpdRel = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax);
primary1IncrSpd = OLMaximizeSpd(cal, primary1IncrSpdRel,'lambda',lambda);
[~,~,primary1IncrSpdPredicted] = OLSpdToSettings(cal,primary1IncrSpd+darkSpd,'lambda', lambda);
primary1IncrSpdPredicted = primary1IncrSpdPredicted-darkSpd;

primary2IncrSpdRel = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax);
primary2IncrSpd = OLMaximizeSpd(cal, primary2IncrSpdRel,'lambda',lambda);
[~,~,primary2IncrSpdPredicted] = OLSpdToSettings(cal,primary2IncrSpd+darkSpd,'lambda', lambda);
primary2IncrSpdPredicted = primary2IncrSpdPredicted-darkSpd;

%% Get set of test lights varying in intensity and primary lights varying
%  in contribution of the two primaries. Scale by the factors defined above
for i = 1:adjustmentLength
    testSpdsNominal(:,i) = testScaleFactor * testScales(i) * testIncrSpd...
        + darkSpd;
    [testSettings(:,i),~,testSpdsPredicted(:,i)] = OLSpdToSettings(cal,...
        testSpdsNominal(:,i), 'lambda', lambda);
    [testStartStops(i,1,:),testStartStops(i,2,:)] = ...
        OLSettingsToStartsStops(cal, testSettings(:,i));
    
    primarySpdsNominal(:,i) = p1ScaleFactor * p1Scales(i) * primary1IncrSpd...
        + p2ScaleFactor * p2Scales(i) * primary2IncrSpd + darkSpd;
    [primarySettings(:,i),~,primarySpdsPredicted(:,i)] = ...
        OLSpdToSettings(cal, primarySpdsNominal(:,i), 'lambda', lambda);
    [primaryStartStops(i,1,:), primaryStartStops(i,2,:)] =...
        OLSettingsToStartsStops(cal, primarySettings(:,i));
end
clear i;

%% Take a look at spectra (optional)
makePlots = false;
if makePlots
    OLPlotSpdCheck(cal.computed.pr650Wls, testSpdsPredicted);
    title('Test');
    
    OLPlotSpdCheck(cal.computed.pr650Wls, primarySpdsPredicted);
    title('Primaries');
end

%% Save results
file = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
    adjustmentLength, p1, p2, test, p1ScaleFactor, p2ScaleFactor, testScaleFactor);
save(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', file));
end