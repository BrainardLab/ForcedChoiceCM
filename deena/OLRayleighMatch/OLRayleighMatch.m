function matches = OLRayleighMatch(varargin)
% Program to run a Rayleigh match experiment on the OneLight
%
% Syntax:
%   OLRayleighMatch
%
% Description:
%    Displays a mixture of two primary iights followed by a test light in
%    alternation.
%
%    Default wavelengths are 580 for test, 670 for p1, and 540 for
%    p2, but these can also be entered by the user as key/value pairs.
%
%    Subjects can use the game pad's directional pad to adjust the
%    intensity of the test light and the ratio of the two primaries to try
%    and get the two fields to match:
%       right moves the primary ratio towards p1 and left towards p2
%       up increases the test intensity and down decreases the test intensity.
%
%    The subject can also use the top 'Y' button to toggle step size, the
%    right 'B' button to record a match, and the bottom 'A' button to quit.
%
%    The routine prompts for subject and session info.
%
% Inputs:
%    None
%
% Outputs:
%    Saves a .mat file titled by subject and session number, which includes
%    a record of subjects' matches and various experimental parameters. As
%    of now, the file is only saved if the subject makes at least one
%    match.
%
% Optional key-value pairs:
%    'p1'        - integer wavelength of the first primary light in nm.
%                  Default is 670.
%    'p2'        - integer wavelength of the second primary light in nm.
%                  Default is 540.
%    'test'      - integer wavelength of the test light in nm. Default is
%                  580.
%    'sInterval' - length of time that each stimulus is displayed for, in
%                   s. Default is 1.5.
%    'isi'       - inter-stimulus interval for white light between primary
%                  and test fields, in s. Default is 0.3.
%    'iti'       - inter-trial interval for white light between test and
%                  primary fields. Default is 1s.
%
%
%    'plotSpds'  - logical indicating whether to make plots of nominal spds.
%                  Default is false.

% History:
%   xx/xx/19  dce       Wrote it.
%   01/15/20  dce, dhb  Add control of white by specifying primaries.
%   01/16/20  dce       Add control of stimulus length and inter-stimulus
%                       intervals

%% Close any stray figures
close all;

%% Parse input
%
% This currently sets the primary and test wavelengths.
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 540, @(x) (isnumeric(x)));
p.addParameter('test', 580, @(x) (isnumeric(x)));
p.addParameter('sInterval', 1.5, @(x) (isnumeric(x)));
p.addParameter('isi', 0.3, @(x) (isnumeric(x)));
p.addParameter('iti', 1, @(x) (isnumeric(x)));
p.addParameter('plotSpds', false, @(x) (islogical(x)));
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
sInterval = p.Results.sInterval;
isi = p.Results.isi;
iti = p.Results.iti;p1 = p.Results.p1;

%% Parameters
%
% Length of scaling factor adjustment arrays
primaries_length = 101;
test_length = 101;

% Scale primary mixture down from max available by
% this factor.
primaryScaleFactor = 1;
if (primaryScaleFactor > 1)
    error('Do not set primaryScaleFactor greater than 1')
end

% Spectrum-generating parameters
fullWidthHalfMax = 20;
lambda = 0.001;

% Scaling factors for intensity of primary and test lights. The two
% primaries are each scaled from 0 to 1, and the arrays are set up so the
% two scaling factors will always add up to 1. We start with primary2
% scaled up to 1 and primary1 scaled down to 0.
p1Scales = linspace(0,1,primaries_length);
p2Scales = 1-p1Scales;
testScales = linspace(0,1,test_length);

%% Set up directory for saving results
subjectID = input('Enter subject ID: ');
sessionNum = input('Enter session number: ');

% Create directory named subjectID for saving data, if it doesn't exist
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Create data file with name subject ID_SessionNumber. Throw error if a
% file already exists with that name.
fileName = [subjectID, '_', num2str(sessionNum), '.mat'];
fileLoc = fullfile(outputDir,fileName);
if (exist(fileLoc, 'file'))
    error('Specified output file %s already exists', fileName);
end


%% Get the OneLight calibration structure
%
% Use it to set some device related parameters.
% numCols is number of mirrors in each OL column
cal = OLGetCalibrationStructure;
[spdLength,settingsLength] = size(cal.computed.pr650M);
numCols = cal.describe.numColMirrors;

%% Initialize arrays for storing precomputed spectra
% These will be cycled through in the adjustments.
% Start-stop arrays hold start positions in the first column and stop
% postions in the second column.
testSpdsNominal = zeros(spdLength,test_length);
testSpdsPredicted = zeros(spdLength,test_length);
testSettings = zeros(settingsLength,test_length);
testStartStops = zeros(test_length,2,numCols);

primarySpdsNominal = zeros(spdLength, primaries_length);
primarySpdsPredicted = zeros(spdLength,primaries_length);
primarySettings = zeros(settingsLength, primaries_length);
primaryStartStops = zeros(primaries_length,2,numCols);

%% Specify and initialize "white" light
%
% This is displayed between the times when subject is comparing
% the primaries and test.
%
% The easiest way to set up something reasonable is to turn
% the OneLight on to half of its max.  We could get fancier
% later and explicitly provide a spectrum.
whitePrimaries = 0.5 * ones(settingsLength, 1);
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
primary1IncrSpdRel = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax);
primary1IncrSpd = OLMaximizeSpd(cal, primary1IncrSpdRel,'lambda',lambda);
primary2IncrSpdRel = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax)/3;
primary2IncrSpd = OLMaximizeSpd(cal, primary2IncrSpdRel,'lambda',lambda);

%% Get set of test lights, varying in intensity
for i = 1:test_length
    testSpdsNominal(:,i) = (testScales(i) * testIncrSpd) + darkSpd;
    [testSettings(:,i),~,testSpdsPredicted(:,i)] = OLSpdToSettings(cal, testSpdsNominal(:,i), 'lambda', lambda);
    [testStartStops(i,1,:),testStartStops(i,2,:)] = OLSettingsToStartsStops(cal, testSettings(:,i));
end

%% Get set of primary lights, varying in relative contribution of two primaries
for i = 1:primaries_length
    primarySpdsNominal(:,i) = (primaryScaleFactor*(p1Scales(i) * primary1IncrSpd) + (p2Scales(i) * primary2IncrSpd)) + darkSpd;
    [primarySettings(:,i),~,primarySpdsPredicted(:,i)] = OLSpdToSettings(cal, primarySpdsNominal(:,i), 'lambda', lambda);
    [primaryStartStops(i,1,:), primaryStartStops(i,2,:)] = OLSettingsToStartsStops(cal, primarySettings(:,i));
end

%% Take a look at spectra (optional)
if p.Results.plotSpds
    figure; clf; hold on; title('Test');
    OLplotSpdCheck(cal.computed.pr650Wls, testSpdsNominal);
    
    figure; clf; title('Primaries');
    OLplotSpdCheck(cal.computed.pr650Wls, primarySpdsNominal);
end

%% Intialize OneLight and button box
ol = OneLight;
gamePad = GamePad();

%% Set up projector
fprintf('\n**** Set up projector ****\n');
annulusFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'), 'projectorSettings','OLAnnulusSettings.mat');
if exist(annulusFile, 'file')
    fprintf('[1]: Use existing annulus settings file\n');
    fprintf('[2]: Reset annulus\n');
    res = GetInput('Select option', 'number', 1);
    if res == 2
        ol.setMirrors(squeeze(primaryStartStops(1,1,:))',...
            squeeze(primaryStartStops(1,2,:))');
        GLW_AnnularStimulusButtonBox();
    end
else
    ol.setMirrors(squeeze(primaryStartStops(1,1,:))',...
        squeeze(primaryStartStops(1,2,:))');ol.setAll(false);
    GLW_AnnularStimulusButtonBox();
end
annulusData = load(annulusFile);
annulusData.win.open;
annulusData.win.draw;
fprintf('\nProjector ready. Starting display loop\n')

%% Display loop
%
% Display parameters
isPrimary = true;           % Are we currently displaying primary or test light?
stepModes = [20 5 1];       % Possible step sizes (relative to adjustment_length)

% Hold information
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array

% Initial position in primary, test, and step size arrays
%
% Start with largest step size
primaryPos = 1;
testPos = 1;
stepModePos = 1;

% Loop through primary and test light until the user presses a key
stillLooping = true;
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
    while(mglGetSecs < nowTime + sInterval)
        key = gamePad.getKeyEvent();
        if (~isempty(key))
            switch(key.charCode)
                case 'GP:Y' % Switch step size mode
                    Snd('Play',sin(0:5000)/50);
                    stepModePos = stepModePos + 1;
                    if stepModePos > length(stepModes)
                        stepModePos = 1;
                    end
                    fprintf('User switched step size to %g \n',...
                        (stepModes(stepModePos)/100.0));
                case 'GP:B' % Subject found a match
                    Snd('Play',sin(0:5000)/50);
                    fprintf('User found match at %g test, %g primary \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    % Save match 
                    matches = [matches;...
                        [testScales(testPos), p1Scales(primaryPos)]];
                    matchPositions = [matchPositions; [testPos, primaryPos]];
                    % Randomize primary and test lights 
                    primaryPos = randi(primaries_length); 
                    testPos = randi(test_length); 
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
    
    % Display "white" light in between iterations. This light is displayed
    % for a short time between primary and test fields (isi) and for a
    % longer time between test and primary fields (iti).
    currTime = mglGetSecs;
    if isPrimary
        delay = isi;
    else
        delay = iti;
    end
    while(mglGetSecs < currTime + delay)
        ol.setMirrors(whiteStarts,whiteStops);
    end
    
    % Switch from primary to test
    isPrimary = ~isPrimary;
end

%% Save matches
if ~isempty(matches)
    save(fileLoc, 'matches', 'matchPositions', 'p1', 'p2', 'test', 'cal',...
        'primarySpdsNominal', 'primarySpdsPredicted', 'testSpdsNominal',...
        'testSpdsPredicted', 'primaryStartStops', 'testStartStops',...
        'whitePrimaries', 'whiteSettings', 'whiteStarts', 'whiteStops',...
        'whiteSpdNominal', 'subjectID', 'sessionNum','annulusData', ...
        'sInterval', 'isi', 'iti');
end

%% Close up
GLW_CloseAnnularStimulus();
ol.setAll(false);

end