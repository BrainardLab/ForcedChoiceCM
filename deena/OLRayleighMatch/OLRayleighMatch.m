function OLRayleighMatch(varargin)
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
%    right 'B' button to record a match, the left 'X' button to change the
%    order that primary and test lights are displayed, and the bottom 'A'
%    button to quit.
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
%                   s. Default is 0.25.
%    'isi'       - inter-stimulus interval for white light between primary
%                  and test fields, in s. Default is 0.25.
%    'iti'       - inter-trial interval for white light between test and
%                  primary fields, in s. Default is 1.
%

% History:
%   xx/xx/19  dce       Wrote it.
%   01/15/20  dce, dhb  Add control of white by specifying primaries.
%   01/16/20  dce       Add control of stimulus length and inter-stimulus
%                       intervals
%   02/4/20   dce       Separated light spectrum computations to a
%                       different file

%% Close any stray figures
close all;

%% Parse input
%
% This currently sets the primary and test wavelengths.
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 540, @(x) (isnumeric(x)));
p.addParameter('test', 580, @(x) (isnumeric(x)));
p.addParameter('sInterval', 0.25, @(x) (isnumeric(x)));
p.addParameter('isi', 0.25, @(x) (isnumeric(x)));
p.addParameter('iti', 1, @(x) (isnumeric(x)));
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
sInterval = p.Results.sInterval;
isi = p.Results.isi;
iti = p.Results.iti;p1 = p.Results.p1;

% Find precomputed spectra, or compute if they do not exist
if (p1 == 670 && p2 == 540 && test == 580)
    fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops', 'OLRayleighMatchFineSpectralSettings.mat');
else
    file = sprintf('OLRayleighMatchSpectralSettings_%g_%g_%g.mat', p1, p2, test);
    fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops', file);
    if ~exist(fName, 'file')
        OLRayleighMatchLightSettings(p1,p2,test);
    end
end

% Load light settings and save needed variables locally
lightSettings = load(fName); 

cal = lightSettings.cal; 
primarySpdsNominal = lightSettings.primarySpdsNominal;
primarySpdsPredicted = lightSettings.primarySpdsPredicted; 
testSpdsNominal = lightSettings.testSpdsNominal;
testSpdsPredicted = lightSettings.testSpdsPredicted; 
primaryStartStops = lightSettings.primaryStartStops; 
testStartStops = lightSettings.testStartStops;
p1Scales = lightSettings.p1Scales; 
testScales = lightSettings.testScales;
whitePrimaries = lightSettings.whitePrimaries; 
whiteSettings = lightSettings.whiteSettings;
whiteStarts = lightSettings.whiteStarts; 
whiteStops = lightSettings.whiteStops; 
whiteSpdNominal = lightSettings.whiteSpdNominal; 
adjustment_length = lightSettings.adjustment_length; 

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
% Possible step sizes (relative to adjustment_length)
stepModes = [floor(adjustment_length/5), floor(adjustment_length/20),...
    floor(adjustment_length/100), floor(adjustment_length/200)];  
lightTimes = [sInterval isi sInterval iti]; % Times for each light settings

% The experiment includes an option to switch the order that primary and
% test lights are displayed. If rev is set to true, lights will be
% displayed in the order specified by lightModeRev instead of the order
% specified by lightMode
rev = false;
lightMode = ['p' 'w' 't' 'w']; % Possible light settings - primary, test, or white
lightModeRev = ['t' 'w' 'p' 'w']; % Switch test and primary order
% lightModeTest = ['t' 'w' 't' 'w']; % Show only test light with white in between

% Data-storing arrays 
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array

% Initial position in primary, test, step size, and light mode arrays
%
% Start with largest step size and with the primary light
primaryPos = 1;
testPos = 1;
stepModePos = 1;
lightModePos = 1;

% Loop through primary and test light until the user presses a key
stillLooping = true;
while(stillLooping)
    nowTime = mglGetSecs;
    
    % Display primary, test, or white light. The white light is displayed
    % for a short time between primary and test lights and a long time
    % between test and primary lights.
        if rev
            lights = lightModeRev;
        else
            lights = lightMode;
        end
    % lights = lightModeTest;
    switch lights(lightModePos)
        case 'p'
            ol.setMirrors(squeeze(primaryStartStops(primaryPos,1,:))',...
                squeeze(primaryStartStops(primaryPos,2,:))');
        case 't'
            ol.setMirrors(squeeze(testStartStops(testPos,1,:))',...
                squeeze(testStartStops(testPos,2,:))');
        case 'w'
            ol.setMirrors(whiteStarts,whiteStops);
    end
    
    % Until time limit runs out, check for user input
    while(mglGetSecs < nowTime + lightTimes(lightModePos))
        key = gamePad.getKeyEvent();
        if (~isempty(key))
            switch(key.charCode)
                case 'GP:Y' % Switch step size mode
                    stepModePos = stepModePos + 1;
                    if stepModePos > length(stepModes)
                        stepModePos = 1;
                    end
                    Snd('Play',sin((0:5000)* stepModePos / 20));
                    fprintf('User switched step size to %g \n',...
                        (stepModes(stepModePos)/ (adjustment_length - 1)));
                    
                case 'GP:B' % Subject found a match
                    Snd('Play',sin(0:5000));
                    fprintf('User found match at %g test, %g primary \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    % Save match
                    matches = [matches;...
                        [testScales(testPos), p1Scales(primaryPos)]];
                    matchPositions = [matchPositions; [testPos, primaryPos]];
                    
                    % Randomize lights and set step size to largest
                    primaryPos = randi(primaries_length);
                    testPos = randi(test_length);
                    stepModePos = 1;
                    
                case 'GP:A' % Quit
                    Snd('Play',sin(0:5000));
                    fprintf('User exited program \n');
                    stillLooping = false;
                    
                case 'GP:X' % Switch order of primary and test lights
                    rev = ~rev;
                    Snd('Play',sin(0:5000));
                    fprintf('User switched order of primary and test lights\n');
                    
                case 'GP:North' % Scale up test intensity
                    testPos = testPos + stepModes(stepModePos);
                    if testPos > test_length
                        Snd('Play',sin(0:5000));
                        fprintf('User reached upper test limit \n');
                        testPos = test_length;
                    end
                    Snd('Play',sin(0:5000)/100);
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case 'GP:South' % Scale down test intensity
                    testPos = testPos - stepModes(stepModePos);
                    if testPos < 1
                        Snd('Play',sin(0:5000));
                        fprintf('User reached lower test limit \n');
                        testPos = 1;
                    end
                    Snd('Play',sin(0:5000)/100);
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case 'GP:East' % Move towards p1
                    primaryPos = primaryPos + stepModes(stepModePos);
                    if primaryPos > adjustment_length
                        Snd('Play',sin(0:5000));
                        fprintf('User reached upper primary limit \n');
                        primaryPos = adjustment_length;
                    end
                    Snd('Play',sin(0:5000)/100);
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case 'GP:West' % Move towards p2
                    primaryPos = primaryPos - stepModes(stepModePos);
                    if primaryPos < 1
                        Snd('Play',sin(0:5000));
                        fprintf('User reached lower primary limit \n');
                        primaryPos = 1;
                    end
                    Snd('Play',sin(0:5000)/100);
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
            end
        end
    end
    
    
    % Switch to the next light to display
    lightModePos = lightModePos + 1;
    if lightModePos > length(light)
        lightModePos = 1;
    end
end

%% Save matches
if ~isempty(matches)
    save(fileLoc, 'matches', 'matchPositions', 'p1', 'p2', 'test', 'cal',...
        'primarySpdsNominal', 'primarySpdsPredicted', 'testSpdsNominal',...
        'testSpdsPredicted', 'primaryStartStops', 'testStartStops',...
        'whitePrimaries', 'whiteSettings', 'whiteStarts', 'whiteStops',...
        'whiteSpdNominal', 'subjectID', 'sessionNum','annulusData', ...
        'sInterval', 'isi', 'iti', 'adjustment_length');
end

%% Close up
GLW_CloseAnnularStimulus();
ol.setAll(false);
end