% Program to run a Rayleigh match experiment on the OneLight

% Description
%    Displays a yellow test light (580nm) followed by a mixture of two
%    primary lights (540nm, 670nm). Subjects can use keypresses to adjust
%    the intensity of the test light and the ratio of the two primaries to
%    try and get the two fields to match: 'a' increases the test intensity, 
%    'b' decreases the test intensity, 'c' moves the primary ratio towards 
%    red, and 'd' moves the primary ratio towards green. The program exits 
%    when the user presses the space bar.

%% Initial parameters
% Base wavelengths
primary1 = 670; % Red
primary2 = 540; % Green
test = 580;

% Scaling factors for intensity of primary and test lights. The two
% primaries are each scaled from 0 to 1, and the arrays are set up so the
% two scaling factors will always add up to 1. We start with primary2
% scaled up to 1 and primary1 scaled down to 0. 
adjustment_length = 21; % Length of scaling factor adjustment arrays
p1Scales = linspace(0,1,adjustment_length); 
p2Scales = linspace(1,0,adjustment_length); 
testScales = linspace(0,1,adjustment_length);

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;
[spdLength,primaryLength,~] =  size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors; % Number of mirrors in each OL column

% Initialize arrays for storing precomputed adjustments. The StartStops 
% arrays have one column each for start and stop. 
testSpds = zeros(spdLength,adjustment_length);
testSettings = zeros(primaryLength,adjustment_length);
testStartStops = zeros(adjustment_length,2,numCols); 

primarySpds = zeros(spdLength,adjustment_length);
primarySettings = zeros(primaryLength,adjustment_length);
primaryStartStops = zeros(adjustment_length,2,numCols); 

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
isPrimary = true; % are we currently displaying primary or test light? 

% Initial position in primary and test arrays 
primaryPos = 1; 
testPos = 1; 

% Enable character listening and turn on OneLight
% ol = OneLight;
mglDisplayCursor(0);
ListenChar(2);
FlushEvents;
fprintf('Starting display loop \n');

% Loop through primary and test light until the user presses a key
while(true)
    % Display primary or test light
%     if isPrimary
%         ol.setMirrors(squeeze(primaryStartStops(primaryPos,1,:))', squeeze(primaryStartStops(primaryPos,2,:))');
%     else
%         ol.setMirrors(squeeze(testStartStops(testPos,1,:))', squeeze(testStartStops(testPos,2,:))'); 
%     end
    
    mglWaitSecs(delaySecs); % Time delay
    isPrimary = ~isPrimary; % Switch from primary to test
    
    % Check for user input
    if CharAvail
        switch(GetChar)
            case ' ' % Exit program
                break;
            case 'a' % Scale up test intensity
                testPos = testPos + 1;
                if testPos > 21
                    testPos = 1;
                end
            case 'b' % Scale down test intensity
                testPos = testPos - 1;
                if testPos < 1 
                    testPos = 1;
                end
            case 'c' % Move towards p1
                primaryPos = primaryPos + 1;
                if primaryPos > 21
                    primaryPos = 21;
                end
            case 'd' % Move towards p2
                primaryPos = primaryPos - 1;
                if primaryPos < 1
                    primaryPos = 1;
                end
        end
        fprintf('User pressed key. Test intensity = %g, red primary = %g \n', testScales(testPos), p1Scales(primaryPos)); 
    end
end

% Clean up once user exits
fprintf('User exited the program \n');
ListenChar(0);
mglDisplayCursor(1);