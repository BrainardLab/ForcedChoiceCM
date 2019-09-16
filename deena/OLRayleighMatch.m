% Program to run a Rayleigh match experiment on the OneLight
% Description
%    Displays a yellow test light (580nm) followed by a mixture of two
%    primary lights (540nm, 670nm). Subjects can use keypresses to adjust
%    the intensity of the test light and the ratio of the two primaries to
%    try and get the two fields to match. The program exits when the user
%    presses the space bar.

%% Initial parameters
% Get the calibration structure
cal = OLGetCalibrationStructure;

% Base wavelengths
test = 580;
primary1 = 670;
primary2 = 540;

% Initial scaling factors for intensity of primary and test lights
p1Scale = 0; %scaling of the red primary, ranging from 0 (green) to 1 (red)
p2Scale = 1 - p1Scale;
testScale = 1;

fullWidthHalfMax = 20;
lambda = 0.001;
[spdLength,primaryLength,~] =  size(cal.computed.pr650M);

% Convert scaled primaries to OL spectra
testSpd = testScale * OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)/3;
testPrimaries = OLSpdToPrimary(cal, testSpd, 'lambda', lambda);
testSettings =  OLPrimaryToSettings(cal, testPrimaries);
[testStart, testStop] = OLSettingsToStartsStops(cal, testSettings);


primary1Spd = OLMakeMonochromaticSpd(cal, primary1, fullWidthHalfMax)/3;
primary2Spd = OLMakeMonochromaticSpd(cal, primary2, fullWidthHalfMax)/3;
totalPrimarySpd = (p1Scale * primary1Spd) + (p2Scale * primary2Spd);

primaryPrimaries = OLSpdToPrimary(cal, totalPrimarySpd, 'lambda', lambda);
primarySettings =  OLPrimaryToSettings(cal, primaryPrimaries);
[primaryStart, primaryStop] = OLSettingsToStartsStops(cal, primarySettings);

%% Display loop
% Display parameters
delaySecs = 2; % time in seconds that a given field is displayed for
isPrimary = true; % logical indicating whether we're currently displaying primary or test light

% Enable character listening and turn on OneLight
mglDisplayCursor(0);
ListenChar(2);
FlushEvents;
ol = OneLight;
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