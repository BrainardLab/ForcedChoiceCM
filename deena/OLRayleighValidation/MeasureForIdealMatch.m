%% Setup
% Load files
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', 'OLRayleighMatchFineSpectralSettings_670_560_600.mat');
lightSettings = load(fName);

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

% Set up initial arrays
ol = OneLight();
[spdLength, ~] = size(lightSettings.primarySpdsPredicted);
primaryData = zeros(spdLength, 21);
testData = zeros(spdLength, 41);

%% Loop
for p = 181:201
    primaryStarts = lightSettings.primaryStartStops(p,1,:);
    primaryStops = lightSettings.primaryStartStops(p,2,:);
    
    % Display primary on OL, measure with radiometer
    ol.setMirrors(squeeze(primaryStarts)', squeeze(primaryStops)');
    pause(0.1);
    primaryMeas = spectroRadiometerOBJ.measure;
    primaryData(:, p - 180) = primaryMeas;
    fprintf('Primary Light %g Complete\n', p-180);
end

for t = 1:41
    testStarts =  lightSettings.testStartStops(t,1,:);
    testStops =  lightSettings.testStartStops(t,2,:);
    
    % Display test on OL, measure with radiometer
    ol.setMirrors(squeeze(testStarts)',...
        squeeze(testStops)');
    pause(0.1);
    testMeas = spectroRadiometerOBJ.measure;
    testData(:,t) = testMeas;
    fprintf('Test Light %g Complete\n', t);
end

% Save data and close devices
fName = '/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/ideal_meas.mat';
save(fName, 'primaryData', 'testData');
spectroRadiometerOBJ.shutDown;
ol.setAll(false);
emailRecipient = 'delul@sas.upenn.edu';
SendEmail(emailRecipient, 'OneLight Measurements Complete', ...
    'Finished!');