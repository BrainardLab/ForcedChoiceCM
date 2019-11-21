function OLRadiometerMatchesPlayback(matches, varargin)

%% Set up initial parameters 
% Base wavelengths - parse input
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 540, @(x) (isnumeric(x)));
p.addParameter('test', 580, @(x) (isnumeric(x)));
p.parse(varargin{:});
p1 = p.Results.p1;
p2 = p.Results.p2; 
test = p.Results.test; 

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
    ol.setMirrors(squeeze(primaryStartStops(i,1,:))',...
            squeeze(primaryStartStops(i,2,:))');
    primaryMeas = spectroRadiometerOBJ.measure;
    measuredPrimarySpds(:,i) = primaryMeas;  
    
    ol.setMirrors(squeeze(testStartStops(i,1,:))',...
            squeeze(testStartStops(i,2,:))');
    testMeas = spectroRadiometerOBJ.measure;
    measuredTestSpds(:,i) = testMeas; 
end
testWls = SToWls(spectroRadiometerOBJ.userS);
save('playbackMatches.mat', 'testWls', 'measuredTestSpds', 'measuredPrimarySpds');
spectroRadiometerOBJ.shutDown;

end 