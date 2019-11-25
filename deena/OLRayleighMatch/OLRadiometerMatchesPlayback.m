function OLRadiometerMatchesPlayback(matchFile)

% Load file and check it contains required variables
load(matchFile);
if (exist('p1', 'var') == 0 || exist('p2', 'var') == 0 ||...
        exist('test', 'var') == 0 || exist('matches', 'var') == 0)
    error('Passed file does not contain required variables');
end
load(matchFile, 'test');  

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;
[spdLength,primaryLength] =  size(cal.computed.pr650M);
[numMatches, ~] = size(matches); % Number of rows in match array 
numCols = cal.describe.numColMirrors; % Number of mirrors in each OL column

% Initialize arrays for storing precomputed adjustments. The StartStops
% arrays have one column each for start and stop.
testSpds = zeros(spdLength,numMatches);
testSettings = zeros(primaryLength,numMatches);
testStartStops = zeros(numMatches,2,numCols);

primarySpds = zeros(spdLength,numMatches);
primarySettings = zeros(primaryLength,numMatches);
primaryStartStops = zeros(numMatches,2,numCols);

% Scale primaries and convert to OL spectra
for i = 1:numMatches
    testSpds(:,i) = matches(i, 1) * OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)/3;
    testSettings(:,i) = OLSpdToSettings(cal, testSpds(:,i), 'lambda', lambda);
    
    [testStart,testStop] = OLSettingsToStartsStops(cal, testSettings(:,i));
    testStartStops(i,1,:) = testStart;
    testStartStops(i,2,:) = testStop;
    
    primary1Spd = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax)/3;
    primary2Spd = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax)/3;
    primarySpds(:,i) = (matches(i,2) * primary1Spd) + ((1 - matches(i, 2)) * primary2Spd);
    primarySettings(:,i) = OLSpdToSettings(cal, primarySpds(:,i), 'lambda', lambda);
    
    [primaryStart,primaryStop] = OLSettingsToStartsStops(cal, primarySettings(:,i));
    primaryStartStops(i,1,:) = primaryStart;
    primaryStartStops(i,2,:) = primaryStop;
end

%% Set up radiometer 
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

%% Loop through matches 
measuredTestSpds = zeros(spdLength, numMatches);
measuredPrimarySpds = zeros(spdLength, numMatches);
ol = OneLight(); 
for i = 1:numMatches
    % Display primary on OL, measure with radiometer
    ol.setMirrors(squeeze(primaryStartStops(i,1,:))',...
            squeeze(primaryStartStops(i,2,:))');
    pause(0.1);
    primaryMeas = spectroRadiometerOBJ.measure;
    measuredPrimarySpds(:,i) = primaryMeas;  
    
    % Display test on OL, measure with radiometer
    ol.setMirrors(squeeze(testStartStops(i,1,:))',...
            squeeze(testStartStops(i,2,:))');
    pause(0.1);
    testMeas = spectroRadiometerOBJ.measure;
    measuredTestSpds(:,i) = testMeas; 
end

% Save data
testWls = SToWls(spectroRadiometerOBJ.userS);
[~, userID] = system('whoami');
userID = strtrim(userID);
fName = fullfile('/Users',userID, 'Documents/MATLAB/projects/Experiments/ForcedChoiceCM/deena/OLRayleighMatch','OLPlaybackMatches.mat');

save(fName, 'testWls', 'measuredTestSpds', 'measuredPrimarySpds');
spectroRadiometerOBJ.shutDown;
disp('Successfully played back matches'); 
end 