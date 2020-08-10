% Script to loop through a series of radiometer measurements to measure
% lights around subjects' match settings. This allows computation of an
% "ideal" match based on measured data.
% Saves a file of measured primary and test spds to a chosen directory.

% History
%    dce    2/xx/20  - Wrote it
%    dce    4/5/20    - Edited for style, added ongoing saving 

%% Parameters - change as needed 
% Light settings to load
file = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', 'OLRayleighMatchFineSpectralSettings_670_560_600.mat');
lightSettings = load(file);

% Indices of primary and test light to measure. These are around the range 
% of subject matches 
pIndices = 181:201;
tIndices = 1:41; 

% Output directory
fName = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'ideal_meas.mat');

% Email recipient - receives message when measurements are complete 
emailRecipient = 'delul@sas.upenn.edu';

%% Setup 
% Set up radiometer
spectroRadiometerOBJ  = PR670dev(...
    'verbosity',        1, ...       % 1 -> minimum verbosity
    'devicePortString', [] ...       % empty -> automatic port detection)
    );

% Set options Options available for PR670:
spectroRadiometerOBJ.setOptions(...
    'verbosity',        1, ...
    'syncMode',         'OFF', ...      % choose from 'OFF', 'AUTO', [20 400];
    'cyclesToAverage',  1, ...          % choose any integer in range [1 99]
    'sensitivityMode',  'STANDARD', ... % choose between 'STANDARD' and 'EXTENDED'.  'STANDARD': (exposure range: 6 - 6,000 msec, 'EXTENDED': exposure range: 6 - 30,000 msec
    'exposureTime',     'ADAPTIVE', ... % choose between 'ADAPTIVE' (for adaptive exposure), or a value in the range [6 6000] for 'STANDARD' sensitivity mode, or a value in the range [6 30000] for the 'EXTENDED' sensitivity mode
    'apertureSize',     '1 DEG' ...   % choose between '1 DEG', '1/2 DEG', '1/4 DEG', '1/8 DEG'
    );

% Set up OneLight
ol = OneLight();

% Set up initial arrays
[spdLength, ~] = size(lightSettings.primarySpdsPredicted);
primaryData = zeros(spdLength, length(pIndices));
testData = zeros(spdLength, length(tIndices));

%% Loop to display and measure lights 
for p = pIndices  % Primaries
    % Starts/stops positions for OneLight
    primaryStarts = lightSettings.primaryStartStops(p,1,:);
    primaryStops = lightSettings.primaryStartStops(p,2,:);
    
    % Display primary on OL, measure with radiometer
    ol.setMirrors(squeeze(primaryStarts)', squeeze(primaryStops)');
    pause(0.1);
    primaryMeas = spectroRadiometerOBJ.measure;
    primaryData(:, p - pIndices(1) + 1) = primaryMeas;
    save(fName, 'primaryData');
    fprintf('Primary Light %g Complete\n', p - pIndices(1) + 1);
end

for t = tIndices  % Test lights 
    % Starts/stops positions for OneLight
    testStarts =  lightSettings.testStartStops(t,1,:);
    testStops =  lightSettings.testStartStops(t,2,:);
    
    % Display test on OL, measure with radiometer
    ol.setMirrors(squeeze(testStarts)',...
        squeeze(testStops)');
    pause(0.1);
    testMeas = spectroRadiometerOBJ.measure;
    testData(:,t - tIndices(1) + 1) = testMeas;
    save(fName, 'primaryData', 'testData');
    fprintf('Test Light %g Complete\n', t - tIndices(1) + 1);
end

% Close devices and email experimenter
spectroRadiometerOBJ.shutDown;
ol.setAll(false);
SendEmail(emailRecipient, 'OneLight Measurements Complete', ...
    'Finished!');