function matches = OLRayleighMatchDichromat(varargin)
% Runs a OneLight Rayleigh match experiment modeled after the Oculus-HMC
% anomaloscope
%
% Syntax:
%   OLRayleighMatchDichromat
%
% Description:
%    Displays various mixtures of two primary lights followed by a test
%    light in alternation. First runs the mixtures with the primary light
%    first, then with the test light first
%
%    Subjects are shown a series of set primary or test lights and use the
%    directional pad to adjust the other light and try to make a match.
%    "Up" increases the test intensity/amount of p2, and "down" decreases
%    the test intensity/amount of p2. For each trial, the subject can then
%    respond either "match" (right 'B' button) or "no match possible"
%    (left 'X' button).
%
%    Default wavelengths are 600 for test, 670 for p1, and 560 for
%    p2, but these can also be entered by the user as key/value pairs.
%
%    The subject can also use the top 'Y' button to toggle step size and
%    the bottom 'A' button to quit.
%
%    The routine prompts for subject and session info.
%
% Inputs:
%    None
%
% Outputs:
%    Saves a .mat file titled by subject and session number, which includes
%    a record of subjects' matches, non-matches, and various experimental
%    parameters.
%
% Optional key-value pairs:
%    'p1'           - integer wavelength of the first primary light in nm.
%                     Default is 670.
%    'p2'           - integer wavelength of the second primary light in nm.
%                     Default is 540.
%    'test'         - integer wavelength of the test light in nm. Default is
%                     580.
%    'sInterval'    - length of time in s that the short stimulus is displayed
%                     for (or both stimuli when white light is used). Default
%                     is 0.25.
%    'lInterval'    - length of time in s that the long stimulus is displayed
%                     for (or the iti white light if it is used). Default is
%                     1.
%    'foveal'       - logical indicating whether we are making foveal matches,
%                     in which case the annulus is not turned on. Default is
%                     true.
%    'setPrimaries' - logical indicating to run a version of the experiment
%                     where the primary lights are set and the user adjusts
%                     the test light, rather than the other way around.
%                     Default is true.

% History:
%   1/22/20  dce   - Wrote program, modified from OLRayleighMatch
%   2/27/20  dce   - Made revisions to program to bring it in line with
%                    OLRayleighMatch.

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
p.addParameter('lInterval', 1, @(x) (isnumeric(x)));
p.addParameter('foveal', true, @(x) (islogical(x)));
p.addParameter('setPrimaries', true, @(x) (islogical(x)));
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
sInterval = p.Results.sInterval;
lInterval = p.Results.lInterval;
foveal = p.Results.foveal;
setPrimaries = p.Results.setPrimaries;

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
fileName = [subjectID, '_dichromat_', num2str(sessionNum), '.mat'];
fileLoc = fullfile(outputDir,fileName);
if (exist(fileLoc, 'file'))
    error('Specified output file %s already exists', fileName);
end

%% Find light settings
% Find precomputed spectra, or compute if they do not exist
file = sprintf('OLRayleighMatchFineSpectralSettings_%g_%g_%g.mat', p1, p2, test);
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', file);
if ~exist(fName, 'file')
    OLRayleighMatchLightSettings(p1,p2,test);
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
adjustment_length = lightSettings.adjustment_length;

%% Intialize OneLight and button box
ol = OneLight;
gamePad = GamePad();

%% Set up projector (if not making foveal matches)
% Set annulusData so the program will save even if annulus is not used
annulusData = 0;
if ~foveal
    fprintf('\n**** Set up projector ****\n');
    annulusFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'projectorSettings','OLAnnulusSettings.mat');
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
% Display parameters
% Possible step sizes (relative to adjustment_length): 20%, 5%, 1%, and
% 0.5% of total range
stepModes = [floor(adjustment_length/5), floor(adjustment_length/20),...
    floor(adjustment_length/100), floor(adjustment_length/200)];
% For the light that can be adjusted, initially set it to 1/3 maximal power
initialAdjustedPos = floor(adjustment_length/3);

% The experiment includes an option to switch the order that primary and
% test lights are displayed. If rev is set to true, lights will be
% displayed in the order specified by lightModeRev instead of the order
% specified by lightMode. In the course of the experiment, this switching
% happens once all the primary lights have been shown.
rev = false;
lightMode = ['p' 't']; % Possible light settings - primary and test
lightModeRev = ['t' 'p']; % Switch test and primary order
lightTimes = [lInterval sInterval]; % Times for each light - first light is shown for longer

% Set indices for values of controlled light. These are tailored for an
% adjustment length of 201 but can be updated as needed.
fixedPositions = [1, 193, 28, 165, 55, 138, 83, 110, 179, 14, 151,...
    41, 124, 69, 96];
fixedPositionPos = 1;

% Data-storing arrays
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array

nonMatches = [];            % Output array with subject mismatches
nonMatchPositions = [];     % Positions of non-matches in adjustment array

% Initial position in light mode, step size, primary, and test arrays
stepModePos = 1;
lightModePos = 1;

if setPrimaries
    primaryPos = fixedPositions(fixedPositionPos);
    testPos = initialAdjustedPos;
else
    testPos = fixedPositions(fixedPositionPos);
    primaryPos = initialAdjustedPos;
end

% Loop control parameters
%
% After a subject makes a match, future matches are blocked until both the
% primary and test lights have been displayed
stillLooping = true;
blockMatches = false;

% Loop through primary and test light until the user presses a key
while(stillLooping)
    nowTime = mglGetSecs;
    
    % Define which light arrays we are selecting from, based on 
    % settings at the current time point  
    if rev
        lights = lightModeRev;
    else
        lights = lightMode;
    end
    
    % Display primary or test light.
    switch lights(lightModePos)
        case 'p'
            ol.setMirrors(squeeze(primaryStartStops(primaryPos,1,:))',...
                squeeze(primaryStartStops(primaryPos,2,:))');
        case 't'
            ol.setMirrors(squeeze(testStartStops(testPos,1,:))',...
                squeeze(testStartStops(testPos,2,:))');
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
                    if ~blockMatches
                        Snd('Play',sin(0:5000));
                        fprintf('User found match at %g test, %g primary \n',...
                            testScales(testPos), p1Scales(primaryPos));
                        matches = [matches; [testScales(testPos),...
                            p1Scales(primaryPos)]];
                        matchPositions = [matchPositions; [testPos,...
                            primaryPos]];
                        % Save data
                        save(fileLoc, 'matches', 'matchPositions','nonMatches', ...
                            'nonMatchPositions', 'fixedPositions', 'p1', 'p2', 'test', 'cal',...
                            'primarySpdsNominal', 'primarySpdsPredicted', 'testSpdsNominal',...
                            'testSpdsPredicted', 'primaryStartStops', 'testStartStops',...
                            'subjectID', 'sessionNum','annulusData','sInterval', 'lInterval',...
                            'foveal', 'setPrimaries');
                        % Move on to next fixed light in fixedPositions
                        stepModePos = 1;
                        fixedPositionPos = fixedPositionPos + 1;
                        if fixedPositionPos <= length(fixedPositions)
                            if setPrimaries
                                primaryPos = fixedPositions(fixedPositionPos);
                                testPos = initialAdjustedPos;
                            else
                                testPos = fixedPositions(fixedPositionPos);
                                primaryPos = initialAdjustedPos;
                            end
                            blockMatches = true;
                        else      % Reached end of fixed positions array
                            if ~rev  % Set up to do reverse matching
                                fixedPositionPos = 1;
                                rev = true;
                            else     % Close program
                                stillLooping = false;
                                fprintf('\nFinished looping through fixed lights\n');
                            end
                        end
                        lightModePos = 0;
                        break;
                    else
                        fprintf('Matching blocked\n');
                    end
                    
                case 'GP:X' % Subject found a non-match
                    if ~ blockMatches
                        Snd('Play',sin(0:5000));
                        fprintf('User found non-match at %g test, %g primary \n',...
                            testScales(testPos), p1Scales(primaryPos));
                        nonMatches = [nonMatches;...
                            [testScales(testPos), p1Scales(primaryPos)]];
                        nonMatchPositions = [nonMatchPositions; [testPos, primaryPos]];
                        % Save data
                        save(fileLoc, 'matches', 'matchPositions','nonMatches', ...
                            'nonMatchPositions', 'fixedPositions', 'p1', 'p2', 'test', 'cal',...
                            'primarySpdsNominal', 'primarySpdsPredicted', 'testSpdsNominal',...
                            'testSpdsPredicted', 'primaryStartStops', 'testStartStops',...
                            'subjectID', 'sessionNum','annulusData','sInterval', 'lInterval',...
                            'foveal', 'setPrimaries');
                        % Move on to next fixed light in fixedPositions
                        stepModePos = 1;
                        fixedPositionPos = fixedPositionPos + 1;
                        if fixedPositionPos <= length(fixedPositions)
                            if setPrimaries
                                primaryPos = fixedPositions(fixedPositionPos);
                                testPos = initialAdjustedPos;
                            else
                                testPos = fixedPositions(fixedPositionPos);
                                primaryPos = initialAdjustedPos;
                            end
                            blockMatches = true;
                        else  % Reached end of fixed positions array
                            if ~rev   % Set up to do reverse matching
                                fixedPositionPos = 1;
                                rev = true;
                            else      % Close program
                                stillLooping = false;
                                fprintf('\nFinished looping through fixed lights\n');
                            end
                        end
                        lightModePos = 0;
                        break;
                    else
                        fprintf('Matching blocked\n');
                    end
                    
                case 'GP:A' % Quit
                    fprintf('User exited program \n');
                    stillLooping = false;
                    Snd('Play',sin(0:5000));
                    break;
                    
                case 'GP:North' % Scale up test intensity or p2 ratio
                    if setPrimaries
                        testPos = testPos + stepModes(stepModePos);
                    else
                        primaryPos = primaryPos + stepModes(stepModePos);
                    end
                    if testPos > adjustment_length || primaryPos > adjustment_length
                        Snd('Play',sin(0:5000));
                        fprintf('User reached upper limit \n');
                        if setPrimaries
                            testPos = adjustment_length;
                        else
                            primaryPos = adjustment_length;
                        end
                    end
                    Snd('Play',sin(0:5000)/100);
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case 'GP:South' % Scale down test intensity or p2 ratio
                    if setPrimaries
                        testPos = testPos - stepModes(stepModePos);
                    else
                        primaryPos = primaryPos - stepModes(stepModePos);
                    end
                    if testPos < 1 || primaryPos < 1
                        Snd('Play',sin(0:5000));
                        fprintf('User reached lower limit \n');
                        if setPrimaries
                            testPos = 1;
                        else
                            primaryPos = 1;
                        end
                    end
                    Snd('Play',sin(0:5000)/100);
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
            end
        end
    end
    
    % Switch to the next light to display
    lightModePos = lightModePos + 1;
    % Check if you have gone through all lights since the last match. If
    % this is the case, turn off match blocking.
    if lightModePos > length(lights)
        lightModePos = 1;
        blockMatches = false;
    end
end

%% Close up
if ~foveal
    GLW_CloseAnnularStimulus();
end
ol.setAll(false);
end