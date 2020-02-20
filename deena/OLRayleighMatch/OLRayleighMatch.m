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
%    'foveal'    - logical indicating whether we are making foveal matches,
%                  in which case the annulus is not turned on. Default is
%                  false.
%    'noWhite'   - logical indicating to run  a version of the experiment
%                  where the lights  are displayed one after the other
%                  without white in between. Default is false.

% History:
%   xx/xx/19  dce       Wrote it.
%   01/15/20  dce, dhb  Add control of white by specifying primaries.
%   01/16/20  dce       Add control of stimulus length and inter-stimulus
%                       intervals
%   02/4/20   dce       Separated light spectrum computations to a
%                       different file
%   02/14/20  dce       Made option to show matches without white light
%% Close any stray figures
close all;

%% Parse input
%
% This currently sets the primary and test wavelengths.
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 560, @(x) (isnumeric(x)));
p.addParameter('test', 600, @(x) (isnumeric(x)));
p.addParameter('sInterval', 0.25, @(x) (isnumeric(x)));
p.addParameter('isi', 0.25, @(x) (isnumeric(x)));
p.addParameter('iti', 1, @(x) (isnumeric(x)));
p.addParameter('foveal', false, @(x) (islogical(x)));
p.addParameter('noWhite', false, @(x) (islogical(x)));
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
sInterval = p.Results.sInterval;
isi = p.Results.isi;
iti = p.Results.iti;p1 = p.Results.p1;
foveal = p.Results.foveal;
noWhite = p.Results.noWhite;

% Find precomputed spectra, or compute if they do not exist
if (p1 == 670 && p2 == 540 && test == 580)
    fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops', 'OLRayleighMatchFineSpectralSettings.mat');
else
    file = sprintf('OLRayleighMatchFineSpectralSettings_%g_%g_%g.mat', p1, p2, test);
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

%% Set up projector (if not making foveal matches)
% Set annulusData so the program will save even if annulus is not used
annulusData = 0;
if ~foveal
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
end

%% Display loop
%
% Display parameters
% Possible step sizes (relative to adjustment_length)
stepModes = [floor(adjustment_length/5), floor(adjustment_length/20),...
    floor(adjustment_length/100), floor(adjustment_length/200)];

% The experiment includes an option to switch the order that primary and
% test lights are displayed. If rev is set to true, lights will be
% displayed in the order specified by lightModeRev instead of the order
% specified by lightMode
rev = false;
if noWhite
    lightMode = ['p', 't'];
    lightModeRev = ['t', 'p'];
    lightTimes = [iti sInterval];
else
    lightMode = ['p' 'w' 't' 'w']; % Possible light settings - primary, test, or white
    lightModeRev = ['t' 'w' 'p' 'w']; % Switch test and primary order
    lightTimes = [sInterval isi sInterval iti]; % Times for each setting
    % lightModeTest = ['t' 'w' 't' 'w']; % Show only test light with white in between
end

% Settings for ideal match. These were derived from the findNominalMatch
% script.
ideal = false;
if foveal
    pIdealIndex = 198;
    tIdealIndex = 12;
else
    pIdealIndex = 199;
    tIdealIndex = 10;
end

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
    if noWhite && (lightModePos == 2)
        Snd('Play',sin(0:5000));
    end
    switch lights(lightModePos)
        case 'p'
            if ideal
                ol.setMirrors(squeeze(primaryStartStops(pIdealIndex,1,:))',...
                    squeeze(primaryStartStops(pIdealIndex,2,:))');
            else
                ol.setMirrors(squeeze(primaryStartStops(primaryPos,1,:))',...
                    squeeze(primaryStartStops(primaryPos,2,:))');
            end
        case 't'
            if ideal
                ol.setMirrors(squeeze(testStartStops(tIdealIndex,1,:))',...
                    squeeze(testStartStops(tIdealIndex,2,:))');
            else
                ol.setMirrors(squeeze(testStartStops(testPos,1,:))',...
                    squeeze(testStartStops(testPos,2,:))');
            end
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
                    for i = 1:stepModePos
                        Snd('Play',sin(0:5000));
                    end
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
                    primaryPos = randi(adjustment_length);
                    testPos = randi(adjustment_length);
                    stepModePos = 1;
                    
                case 'GP:A' % Quit
                    Snd('Play',sin(0:5000));
                    fprintf('User exited program \n');
                    stillLooping = false;
                    
                    % Switch order of primary and test lights. One beep means
                    % primary is now first, two beeps means test is first.
                case 'GP:X'
                    rev = ~rev;
                    Snd('Play',sin(0:5000)/ 20);
                    if rev
                        Snd('Play',sin(0:5000)/ 20);
                    end
                    fprintf('User switched order of primary and test lights\n');
                    
                case'GP:Back' %Switch to showing ideal match
                    ideal = ~ ideal;
                    Snd('Play',sin(0:5000));
                    
                case 'GP:North' % Scale up test intensity
                    testPos = testPos + stepModes(stepModePos);
                    if testPos > adjustment_length
                        Snd('Play',sin(0:5000));
                        fprintf('User reached upper test limit \n');
                        testPos = adjustment_length;
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
    
    % In mode without white light, play sound before switching
    if noWhite && (lightModePos == 2) 
        Snd('Play',sin(0:5000));
    end
    
    % Switch to the next light to display
    lightModePos = lightModePos + 1;
    if lightModePos > length(lights)
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
        'sInterval', 'isi', 'iti', 'adjustment_length', 'foveal');
end

%% Close up
if ~foveal
    GLW_CloseAnnularStimulus();
end
ol.setAll(false);
end