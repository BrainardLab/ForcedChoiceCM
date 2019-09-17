% Program to run a Rayleigh match experiment on the OneLight

% Description
%    Displays a yellow test light (580nm) followed by a mixture of two
%    primary lights (540nm, 670nm). Subjects can use keypresses to adjust
%    the intensity of the test light and the ratio of the two primaries to
%    try and get the two fields to match. The program exits when the user
%    presses the space bar.

%% Initial parameters
% Base wavelengths
test = 580;
primary1 = 670;
primary2 = 540;

adjustment_length = 21; % length of scaling factor adjustment arrays

% Scaling factors for intensity of primary and test lights
p1Scales = linspace(0,1,adjustment_length); %scaling of the red primary, ranging from 0 (green) to 1 (red)
p2Scales = linspace(1,0,adjustment_length);
testScales = linspace(0,1,adjustment_length);

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;
[spdLength,primaryLength,~] =  size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors; % number of mirrors in each OL column

% Initialize arrays for storing precomputed adjustments 
testSpds = zeros(spdLength,adjustment_length);
testSettings = zeros(primaryLength,adjustment_length);
testStartStops = zeros(adjustment_length,2,numCols); % one column each for start and stop

primarySpds = zeros(spdLength,adjustment_length);
primarySettings = zeros(primaryLength,adjustment_length);
primaryStartStops = zeros(adjustment_length,2,numCols); % one column each for start and stop

% Scale primaries and convert to OL spectra
for i = 1:adjustment_length
    testSpds(:,i) = testScales(i) * OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)/3;
    testSettings(:,i) = OLSpdToSettings(cal, testSpds(:,i), 'lambda', lambda);
    
    [testStart,testStop] = OLSettingsToStartsStops(cal, testSettings(:,i));
    testStartStops(i,1,:) = testStart;
    testStartStops(i,2,:) = testStop;
     
    primary1Spd = OLMakeMonochromaticSpd(cal, primary1, fullWidthHalfMax)/3;
    primary2Spd = OLMakeMonochromaticSpd(cal, primary2, fullWidthHalfMax)/3;
    primarySpds(:,i) = (p1Scales(i) * primary1Spd) + (p2Scales(i) * primary2Spd);
    primarySettings(:,i) = OLSpdToSettings(cal, primarySpds(:,i), 'lambda', lambda);
    
    [primaryStart,primaryStop] = OLSettingsToStartsStops(cal, primarySettings(:,i));
    primaryStartStops(i,1,:) = primaryStart;
    primaryStartStops(i,2,:) = primaryStop;
end

%% Display loop
% Display parameters
delaySecs = 2; % time in seconds that a given field is displayed for
isPrimary = true; % logical indicating whether we're currently displaying primary or test light

% Initial position in primary and test arrays 
primary_pos = 1; 
test_pos = 1; 

% Enable character listening and turn on OneLight
ol = OneLight;
mglDisplayCursor(0);
ListenChar(2);
FlushEvents;
fprintf('Starting display loop \n');

% Loop through primary and test light until the user presses a key
while(true)
    % display primary or test light
    if isPrimary
        ol.setMirrors(primaryStart, primaryStop);
    else
        ol.setMirrors(testStart, testStop);
    end
    mglWaitSecs(delaySecs); % time delay
    isPrimary = ~isPrimary; % switch from primary to test
    if CharAvail
        switch(GetChar)
            case ' '
                break;
            case 'a' % scale up test intensity
                testScale = testScale + 0.1;
                if testScale > 1
                    testScale = 1;
                end
            case 'b' % scale down test intensity
                testScale = testScale - 0.1;
                if testScale < 0
                    test1Scale = 0;
                end
            case 'c' % move towards red
                p1Scale = p1Scale + 0.1;
                if p1Scale > 1
                    p1Scale = 1;
                end
            case 'd' % move towards green
                p1Scale = p1Scale - 0.1;
                if p1Scale < 0
                    p1Scale = 0;
                end
        end
    end
end

% Clean up once user exits
fprintf('User exited the program \n');
ListenChar(0);
mglDisplayCursor(1);