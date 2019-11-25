
% Program to run a Rayleigh match experiment on the OneLight
function matches = OLRayleighMatch(varargin)
% Syntax:
%   OLRayleighMatch(varargin)
%
% Description
%    Displays a yellow test light followed by a mixture of two primary
%    lights. Default wavelengths are 580 for test, 670 for p1, and 540 for
%    p2, but these can also be entered by the user/ Subjects can use the
%    game pad's directional pad to adjust the intensity of the test light
%    and the ratio of the two primaries to try and get the two fields to
%    match: right moves the primary ratio towards p1 and left towards
%    p2, up increases the test intensity and down decreases the test
%    intensity. The subject can also use the top 'Y' button to toggle step
%    size, the right 'B' button to record a match, and the bottom 'A'
%    button to quit.
%
% Inputs:
%    none
%
% Outputs:
%    matches - 2x2 integer array containing test and primary values for
%              subject's matches (first column is test, second is primary)
%
% Optional key-value pairs:
%    'p1'    - integer wavelength of the first primary light in nm. Default
%              is 670. 
%    'p2'    - integer wavelength of the second primary light in nm.
%              Default is 540.
%    'test'  - integer wavelength of the test light in nm. Default is 580.



%% Initial parameters
% Base wavelengths - parse input
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 540, @(x) (isnumeric(x)));
p.addParameter('test', 580, @(x) (isnumeric(x)));
p.parse(varargin{:});
p1 = p.Results.p1;
p2 = p.Results.p2; 
test = p.Results.test; 

% Scaling factors for intensity of primary and test lights. The two
% primaries are each scaled from 0 to 1, and the arrays are set up so the
% two scaling factors will always add up to 1. We start with primary2
% scaled up to 1 and primary1 scaled down to 0.
adjustment_length = 101; % Length of scaling factor adjustment arrays
p1Scales = linspace(0,1,adjustment_length);
p2Scales = linspace(1,0,adjustment_length);
testScales = linspace(0,1,adjustment_length);

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;
[spdLength,primaryLength] =  size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors; % Number of mirrors in each OL column

% Initialize arrays for storing precomputed adjustments. The StartStops
% arrays have one column each for start and stop.
testSpdsNominal = zeros(spdLength,adjustment_length);
testSpdsPredicted = zeros(spdLength,adjustment_length);
testSettings = zeros(primaryLength,adjustment_length);
testStartStops = zeros(adjustment_length,2,numCols);

primarySpdsNominal = zeros(spdLength,adjustment_length);
primarySpdsPredicted = zeros(spdLength,adjustment_length);
primarySettings = zeros(primaryLength,adjustment_length);
primaryStartStops = zeros(adjustment_length,2,numCols);

% Scale primaries and convert to OL spectra
for i = 1:adjustment_length
    testSpdsNominal(:,i) = testScales(i) * OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)/3;
    [testSettings(:,i),~,testSpdsPredicted(:,i)] = OLSpdToSettings(cal, testSpdsNominal(:,i), 'lambda', lambda);
    
    % Check whether any of the settings are 1. If so, we can't actually
    % produce the desired spectrum.
    if (any(testSettings(:,i) == 1))
        fprintf('Maxed out on test settings.\n');
    end
    
    [testStart,testStop] = OLSettingsToStartsStops(cal, testSettings(:,i));
    testStartStops(i,1,:) = testStart;
    testStartStops(i,2,:) = testStop;
    
    primary1Spd = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax)/3;
    primary2Spd = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax)/3;
    primarySpdsNominal(:,i) = (p1Scales(i) * primary1Spd) + (p2Scales(i) * primary2Spd);
    [primarySettings(:,i),~,primarySpdsPredicted(:,i)] = OLSpdToSettings(cal, primarySpdsNominal(:,i), 'lambda', lambda);
    
    if (any(primarySettings(:,i) == 1))
        fprintf('Maxed out on primary settings.\n');
    end
    
    [primaryStart,primaryStop] = OLSettingsToStartsStops(cal, primarySettings(:,i));
    primaryStartStops(i,1,:) = primaryStart;
    primaryStartStops(i,2,:) = primaryStop;
end

%% Take a look at spectra (optional)
makeFigs = false;
if makeFigs
    figure; clf; hold on
    OLplotSpdCheck(testSpdsNominal,cal);
    
    figure; clf;
    OLplotSpdCheck(primarySpdsNominal,cal);
end

%% Display loop
% Display parameters
delaySecs = 2; % time in seconds that a given field is displayed for
isPrimary = true; % are we currently displaying primary or test light?
stepModes = [20 5 1]; % Possible step sizes (relative to adjustment_length)
matches = []; % Output array with subject matches

% Initial position in primary, test, and step size arrays
primaryPos = 1;
testPos = 1;
stepModePos = 1; % Start with largest step size

% Intialize OneLight and button box
ol = OneLight;
gamePad = GamePad();
fprintf('Starting display loop \n');
stillLooping = true;

% Loop through primary and test light until the user presses a key
while(stillLooping)
    nowTime = mglGetSecs;
    % Display primary or test light
    if isPrimary
        Snd('Play',sin(0:5000));
        ol.setMirrors(squeeze(primaryStartStops(primaryPos,1,:))',...
            squeeze(primaryStartStops(primaryPos,2,:))');
    else
        Snd('Play',sin((0:5000)/100));
        ol.setMirrors(squeeze(testStartStops(testPos,1,:))',...
            squeeze(testStartStops(testPos,2,:))');
    end
    
    % Until time limit runs out, check for user input
    while(mglGetSecs < nowTime + delaySecs)
        key = gamePad.getKeyEvent();
        if (~isempty(key))
            switch(key.charCode)
                case 'GP:Y' % Switch step size mode
                    stepModePos = stepModePos + 1;
                    if stepModePos > length(stepModes)
                        stepModePos = 1;
                    end
                    fprintf('User switched step size to %g \n',...
                        (stepModes(stepModePos)/100.0));
                case 'GP:B' % Subject found a match
                    fprintf('User found match at %g test, %g primary \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    matches = [matches;...
                        [testScales(testPos), p1Scales(primaryPos)]];
                    matchPositions = [matchPositions; [testPos, primaryPos]]; 
                case 'GP:A' % Quit
                    fprintf('User exited program \n');
                    stillLooping = false; 
                case 'GP:North' % Scale up test intensity
                    testPos = testPos + stepModes(stepModePos);
                    if testPos > adjustment_length
                        testPos = adjustment_length;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                case 'GP:South' % Scale down test intensity
                    testPos = testPos - stepModes(stepModePos);
                    if testPos < 1
                        testPos = 1;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                case 'GP:East' % Move towards p1
                    primaryPos = primaryPos + stepModes(stepModePos);
                    if primaryPos > adjustment_length
                        primaryPos = adjustment_length;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                case 'GP:West' % Move towards p2
                    primaryPos = primaryPos - stepModes(stepModePos);
                    if primaryPos < 1
                        primaryPos = 1;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
            end
        end
    end
    isPrimary = ~isPrimary; % Switch from primary to test
end

% Save matches 
[~, userID] = system('whoami');
userID = strtrim(userID);
fName = fullfile('/Users',userID, 'Documents/MATLAB/projects/Experiments/ForcedChoiceCM/deena/OLRayleighMatch','OLSampleMatches.mat');
save(fName, 'matches', 'matchPositions', 'p1', 'p2', 'test', 'cal', 'primarySpdsNominal', 'primarySpdsPredicted', 'testSpdsNominal', 'testSpdsPredicted', 'primaryStartStops', 'testStartStops');
end