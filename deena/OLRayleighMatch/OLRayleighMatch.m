
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
% Set up director to save results into
subjectID = input('Enter subject ID: ');
sessionNum = input('Enter session number: ');

% Create directory named SubjectID for saving data, if it doesn't exist
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Create data file with name Subject ID_SessionNumber. Throw error if a
% file already exists with that name
fileName = [subjectID, '_', num2str(sessionNum), '.mat'];
fileLoc = fullfile(outputDir,fileName);
if (exist(fileLoc, 'file'))
    error('Specified output file %s already exists', fileName);
end

% Base wavelengths - parse input
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 540, @(x) (isnumeric(x)));
p.addParameter('test', 580, @(x) (isnumeric(x)));
p.parse(varargin{:});
p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;

% Length of scaling factor adjustment arrays
primaries_length = 101;
test_length = 101;

% Scaling factors for intensity of primary and test lights. The two
% primaries are each scaled from 0 to 1, and the arrays are set up so the
% two scaling factors will always add up to 1. We start with primary2
% scaled up to 1 and primary1 scaled down to 0.
p1Scales = linspace(0,1,primaries_length);
p2Scales = linspace(1,0,primaries_length);
testScales = linspace(0,1,test_length);

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;
[spdLength,settingsLength] =  size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors; % Number of mirrors in each OL column

% Initialize arrays for storing precomputed adjustments. The StartStops
% arrays have one column each for start and stop.
testSpdsNominal = zeros(spdLength,test_length);
testSpdsPredicted = zeros(spdLength,test_length);
testSettings = zeros(settingsLength,test_length);

primarySpdsNominal = zeros(spdLength, primaries_length);
primarySpdsPredicted = zeros(spdLength,primaries_length);
primarySettings = zeros(settingsLength, primaries_length);

% Scale test and primary lights, clip if needed, convert to OL spectra
for i = 1:test_length
    testSpdsNominal(:,i) = testScales(i) * OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)/3;
    [testSettings(:,i),~,testSpdsPredicted(:,i)] = OLSpdToSettings(cal, testSpdsNominal(:,i), 'lambda', lambda);
    
    % Check whether any of the settings are 1. If so, we can't actually
    % produce the desired spectrum.
    if testSettings(:,i) == 1
        testSettings = testSettings(:, 1:(i-1));
        testSpdsPredicted = testSpdsPredicted(:, 1:(i-1));
        test_length = i - 1;
        fprintf('Maxed out on test settings at position %g.\n', i);
        break;
    end
end

testStartStops = zeros(test_length,2,numCols);
for i = 1: test_length
    [testStart,testStop] = OLSettingsToStartsStops(cal, testSettings(:,i));
    testStartStops(i,1,:) = testStart;
    testStartStops(i,2,:) = testStop;
end

for i = 1:primaries_length
    primary1Spd = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax)/3;
    primary2Spd = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax)/3;
    primarySpdsNominal(:,i) = (p1Scales(i) * primary1Spd) + (p2Scales(i) * primary2Spd);
    [primarySettings(:,i),~,primarySpdsPredicted(:,i)] = OLSpdToSettings(cal, primarySpdsNominal(:,i), 'lambda', lambda);
    if primarySettings(:,i) == 1
        primarySettings = primarySettings(:, 1:(i-1));
        primarySpdsPredicted = primarySpdsPredicted(:, 1:(i-1));
        primaries_length = i - 1;
        fprintf('Maxed out on primary settings at position %g.\n', i);
        break;
    end
end
primaryStartStops = zeros(primaries_length,2,numCols);
for i = 1:primaries_length
    [primaryStart,primaryStop] = OLSettingsToStartsStops(cal, primarySettings(:,i));
    primaryStartStops(i,1,:) = primaryStart;
    primaryStartStops(i,2,:) = primaryStop;
end
%% Take a look at spectra (optional)
makeFigs = false;
if makeFigs
    figure; clf; hold on
    OLplotSpdCheck(cal.computed.pr650Wls, testSpdsNominal);
    
    figure; clf;
    OLplotSpdCheck(cal.computed.pr650Wls, primarySpdsNominal);
end

%% Set up projector
fprintf('**** Set up projector ****\n'); 
annulusFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'OLAnnulusSettings.mat'); 
if exist(annulusFile, 'file')
    fprintf('\n'); 
    fprintf(['[1]: Use existing annulus settings file ', annulusFile, '\n']);
    fprintf('[2]: Reset annulus\n');
    res = GetInput('Select option:', 'number', 1);
    if res == 1 
        annulusData = load(annulusFile); 
        annulusData.win.open; 
        annulusData.win.draw; 
    else 
        GLW_AnnularStimulusButtonBox(); 
    end 
else 
    GLW_AnnularStimulusButtonBox(); 
end 
fprintf('\nProjector ready. Starting display loop\n') 

%% Display loop
% Display parameters
delaySecs = 2; % time in seconds that a given field is displayed for
isPrimary = true; % are we currently displaying primary or test light?
whiteStopIndex = 700; % start stops index for brightness of white light 
stepModes = [20 5 1]; % Possible step sizes (relative to adjustment_length)
matches = []; % Output array with subject matches
matchPositions = []; % Positions o matches in the adjustment array

% Initial position in primary, test, and step size arrays
primaryPos = 1;
testPos = 1;
stepModePos = 1; % Start with largest step size

% Intialize OneLight and button box
ol = OneLight;
gamePad = GamePad();
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
                    if testPos > test_length
                        Snd('Play',sin(0:5000)/50);
                        fprintf('User reached upper test limit \n');
                        testPos = test_length;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                case 'GP:South' % Scale down test intensity
                    testPos = testPos - stepModes(stepModePos);
                    if testPos < 1
                        Snd('Play',sin(0:5000)/50);
                        fprintf('User reached lower test limit \n');
                        testPos = 1;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                case 'GP:East' % Move towards p1
                    primaryPos = primaryPos + stepModes(stepModePos);
                    if primaryPos > primaries_length
                        Snd('Play',sin(0:5000)/50);
                        fprintf('User reached upper primary limit \n');
                        primaryPos = primaries_length;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                case 'GP:West' % Move towards p2
                    primaryPos = primaryPos - stepModes(stepModePos);
                    if primaryPos < 1
                        Snd('Play',sin(0:5000)/50);
                        fprintf('User reached lower primary limit \n');
                        primaryPos = 1;
                    end
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
            end
        end
    end
    % Display white light for one second on every other iteration
        if ~isPrimary 
            nowTime = mglGetSecs;
            while(mglGetSecs < nowTime + 1)
                ol.setMirrors(squeeze(whiteStopIndex * ones(1, numCols))',...
                squeeze(whiteStopIndex * ones(1, numCols))');
            end 
        end
    isPrimary = ~isPrimary; % Switch from primary to test
end

GLW_CloseAnnularStimulus();
ol.setAll(false); 
% Save matches
if ~isempty(matches)
    save(fileLoc, 'matches', 'matchPositions', 'p1', 'p2', 'test', 'cal',...
    'primarySpdsNominal', 'primarySpdsPredicted', 'testSpdsNominal',...
    'testSpdsPredicted', 'primaryStartStops', 'testStartStops',...
    'subjectID', 'sessionNum');
end 
end