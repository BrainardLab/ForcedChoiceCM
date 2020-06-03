function OLRayleighMatch(varargin)
% Program to run a Rayleigh match experiment on the OneLight, or a
% simulated Rayleigh match experiment
%
% Syntax:
%   OLRayleighMatch
%
% Description:
%    Displays a mixture of two primary iights followed by a test light in
%    alternation.
%
%    Default wavelengths are 600 for test, 670 for p1, and 560 for
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
%    order that primary and test lights are displayed, the bottom 'A'
%    button to quit, and the 'Back' button to display an ideal match
%    (based on measured data).
%
%    In addition to the experimental program, two simulation options are
%    available. The first option lets the experimenter use keypresses to
%    recreate subjects' button presses. The second option creates a
%    simulated observer based on the Asano model and searches for the
%    best match for this observer.
%
%    The routine prompts for subject and session info.
%
% Inputs:
%    None
%
% Outputs:
%    Saves a .mat file titled by subject and session number, which includes
%    a record of subjects' matches and various experimental parameters. The
%    file is saved after each match and is only saved if the subject makes
%    at least one match.
%
% Optional key-value pairs:
%    'p1'        - integer wavelength of the first primary light in nm.
%                  Default is 670.
%    'p2'        - integer wavelength of the second primary light in nm.
%                  Default is 540.
%    'test'      - integer wavelength of the test light in nm. Default is
%                  580.
%    'sInterval' - length of time in s that the short stimulus is displayed
%                  for (or both stimuli when white light is used). Default
%                  is 0.25.
%    'lInterval' - length of time in s that the long stimulus is displayed
%                  for (or the iti white light if it is used). Default is
%                  1.
%    'foveal'    - logical indicating whether we are making foveal matches,
%                  in which case the annulus is not turned on. Default is
%                  true.
%    'white'     - logical indicating to run  a version of the experiment
%                  where the lights  are displayed with white in between.
%                  Default is false.
%    'simulate'  - logical. Set to true to run in simulated mode. Default
%                  is true.
%    'silent'    - logical. Set to true to turn off feedback beeps. Default
%                  is false.

% History:
%   xx/xx/19  dce       Wrote it.
%   01/15/20  dce, dhb  Add control of white by specifying primaries.
%   01/16/20  dce       Add control of stimulus length and inter-stimulus
%                       intervals
%   02/4/20   dce       Separated light spectrum computations to a
%                       different file
%   02/14/20  dce       Made option to show matches without white light
%   02/26/19  dce       Changed default settings to foveal and no white
%                       light
%   03/29/19  dce       Edited for style
%   06/01/20  dce       Added simulation option
%   06/02/20  dce       Added simulated observer logic
%   06/03/20  dce       Added tracking of subject settings over time

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
p.addParameter('white', false, @(x) (islogical(x)));
p.addParameter('simulate', true, @(x) (islogical(x)));
p.addParameter('silent', false, @(x) (islogical(x)));

p.parse(varargin{:});
p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
sInterval = p.Results.sInterval;
lInterval = p.Results.lInterval;
foveal = p.Results.foveal;
white = p.Results.white;
simulate = p.Results.simulate;
silent = p.Results.silent;

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

%% Initialize simulated observer if desired
% If the experiment is run in simulation mode without a simulated observer,
% the user can control the stimuli with keypresses
sim_observer = false;
observer = struct();
if simulate
    fprintf('Use simulated observer?\n')
    fprintf('[1]: Yes\n');
    fprintf('[2]: No\n');
    res = GetInput('Select option', 'number', 1);
    if res == 1
        sim_observer = true;
        nMatches = GetInput('Enter number of simulated matches', 'number', 1);
        params = GetInput('Enter optional observer params vector, or press Enter to continue', 'number', -1);
        if isempty(params)
            observer = genRayleighObserver(foveal);
        elseif length(params) == 9
            observer = genRayleighObserver(foveal, 'coneVec', params);
        else
            error('Observer parameters must be listed as a 9-element vector');
        end
    end
end

%% Set up key interpretations
if simulate
    keyCodes.switchStepSize = 's';
    keyCodes.switchLightOrder = 'o';
    keyCodes.foundMatch = ' ';
    keyCodes.quit = 'q';
    keyCodes.idealMatch = 'i';
    keyCodes.increaseIntensity = 'u';
    keyCodes.decreaseIntensity = 'd';
    keyCodes.increaseP1 = 'r';
    keyCodes.decreaseP1 = 'l';
else
    keyCodes.switchStepSize = 'GP:Y';
    keyCodes.switchLightOrder = 'GP:X';
    keyCodes.foundMatch = 'GP:B';
    keyCodes.quit = 'GP:A';
    keyCodes.idealMatch = 'GP:Back';
    keyCodes.increaseIntensity = 'GP:North';
    keyCodes.decreaseIntensity = 'GP:South';
    keyCodes.increaseP1 = 'GP:East';
    keyCodes.decreaseP1 = 'GP:West';
end

%% Find light settings
% Find precomputed spectra, or compute if they do not exist
file = sprintf('OLRayleighMatchFineSpectralSettings_%g_%g_%g.mat', p1, p2, test);
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', file);
if ~exist(fName, 'file')
    OLRayleighMatchLightSettings(p1,p2,test);
end

% Load light settings and save some variables locally
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

%% Intialize OneLight and button box
ol = OneLight('simulate', simulate, 'plotWhenSimulating', true);
if simulate && ~sim_observer
    ListenChar(2);
    FlushEvents;
elseif ~simulate
    gamePad = GamePad();
end

%% Set up projector (if not making foveal matches)
annulusData = struct(); % Empty placeholder variable
% Check if an annulus file exists. If it does, can choose whether to use
% the existing file or reset the annulus.
if ~simulate && ~foveal
    fprintf('\n**** Set up projector ****\n');
    annulusFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'projectorSettings','OLAnnulusSettings.mat');
    if exist(annulusFile, 'file')
        fprintf('[1]: Use existing annulus settings file\n');
        fprintf('[2]: Reset annulus\n');
        res = GetInput('Select option', 'number', 1);
        % Run program to reset the annulus
        if res == 2
            ol.setMirrors(squeeze(primaryStartStops(1,1,:))',...
                squeeze(primaryStartStops(1,2,:))');
            GLW_AnnularStimulusButtonBox();
        end
        % Run program to reset the annulus if its file does not exist
    else
        ol.setMirrors(squeeze(primaryStartStops(1,1,:))',...
            squeeze(primaryStartStops(1,2,:))');ol.setAll(false);
        GLW_AnnularStimulusButtonBox();
    end
    % Load and display annulus file
    annulusData = load(annulusFile);
    annulusData.win.open;
    annulusData.win.draw;
    fprintf('\nProjector ready. Starting display loop\n')
end

%% Display parameters
% Possible step sizes (relative to adjustment_length)
stepModes = [floor(adjustment_length/5), floor(adjustment_length/20),...
    floor(adjustment_length/100), floor(adjustment_length/200)];

% The experiment includes an option to switch the order that primary and
% test lights are displayed. If rev is set to true, lights will be
% displayed in the order specified by lightModeRev instead of the order
% specified by lightMode. When the white light is used, it is displayed for
% a short time between primary and test lights and a long time between test
% and primary lights.
if white
    lightMode = ['p' 'w' 't' 'w']; % Possible light settings - primary, test, or white
    lightModeRev = ['t' 'w' 'p' 'w']; % Switch test and primary order
    lightTimes = [sInterval sInterval sInterval lInterval]; % Times for each setting
else
    lightMode = ['p', 't'];     % Possible light settings - primary, test
    lightModeRev = ['t', 'p'];  % Switch primary and test lights
    lightTimes = [lInterval sInterval]; % Times for each setting
end

% Settings for ideal match. These were derived from the findNominalMatch
% script. Should probably update in the future.
if foveal
    pIdealIndex = 197;
    tIdealIndex = 14;
else
    pIdealIndex = 198;
    tIdealIndex = 13;
end

%% Setup for display loop
% Data-storing arrays
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array

% Initial settings
primaryPos = 1;      % Start with first position in primary array
testPos = 1;         % Start with first position in test array
stepModePos = 1;     % Start with the largest step size
lightModePos = 1;    % Start with the first light in the lights array
rev = false;         % Start with forward, not reverse order
ideal = false;       % Do not start with the ideal match
stillLooping = true; % Start looping

% Store positions of subjects' light settings
subjectSettings = [testScales(testPos), p1Scales(primaryPos)];

% Settings for simulated observer
if sim_observer
    primary_set = false;  % Primary ratio has not yet been adjusted
    test_set = false;     % Test intensity has not yet been adjusted
    
    p1_up_prev = true;    % Indicate that the program directed p1 and t to
    t_up_prev = true;     % increase on prior trials.
end

%% Display loop
% Loop through primary and test light until the user presses a key or the
% simulated observer decides to simulate a keypress
while(stillLooping)
    % Determine which lights we are using
    if rev
        lights = lightModeRev;
    else
        lights = lightMode;
    end
    if ideal
        pI = pIdealIndex;
        tI = tIdealIndex;
    else
        pI = primaryPos;
        tI = testPos;
    end
    
    % Display primary, test, or white light.
    nowTime = mglGetSecs;
    switch lights(lightModePos)
        case 'p'
            ol.setMirrors(squeeze(primaryStartStops(pI,1,:))',...
                squeeze(primaryStartStops(pI,2,:))');
        case 't'
            ol.setMirrors(squeeze(testStartStops(tI,1,:))',...
                squeeze(testStartStops(tI,2,:))');
        case 'w'
            ol.setMirrors(whiteStarts,whiteStops);
    end
    
    % Until time limit runs out, check for user input
    while(mglGetSecs < nowTime + lightTimes(lightModePos))
        % Simulated observer experiment - determine appropriate action with
        % calculations. First adjust the primary ratio, then adjust the
        % test intensity.
        if sim_observer
            if ideal
                pI = pIdealIndex;
                tI = tIdealIndex;
            else
                pI = primaryPos;
                tI = testPos;
            end
            [p1_up, t_up] = observerRayleighDecision(observer,...
                primarySpdsPredicted(:, pI), testSpdsPredicted(:, tI),...
                'noisy', true);
            
            if ~primary_set  % Adjusting primary ratio
                if p1_up == p1_up_prev % Continue in existing direction
                    if p1_up
                        key.charCode = 'r';
                    else
                        key.charCode = 'l';
                    end
                else        % Reversal
                    key.charCode = 's'; % Lower step size
                    % If smallest step size has already been reached,
                    % you're done with the primary.
                    if stepModePos == length(stepModes)
                        primary_set = true;
                    end
                end
                p1_up_prev = p1_up;
                
            elseif ~test_set % Adjusting test intensity
                if t_up == t_up_prev % Continue in existing direction
                    if t_up
                        key.charCode = 'u';
                    else
                        key.charCode = 'd';
                    end
                else        % Reversal
                    % If smallest step size has already been reached,
                    % you're done with the test light and can record a match.
                    if stepModePos == length(stepModes)
                        test_set = true;
                        key.charCode = ' ';
                    else
                        key.charCode = 's'; % Lower step size
                    end
                end
                t_up_prev = t_up;
                
            else
                % Reach this block once a match has been completed. Quit if
                % the desired number of matches has been made, and reset
                % otherwise.
                [row,~] = size(matches);
                if row == nMatches
                    key.charCode = 'q';
                else
                    p1_up_prev = true;
                    t_up_prev = true;
                    primary_set = false;
                    test_set = false;
                    continue;
                end
            end
            
        elseif simulate      % Simulation with user keypresses
            if CharAvail
                key.charCode = GetChar(true,false);
            else
                key = [];
            end
            
        else                 % Live experiment
            key.charCode = gamePad.getKeyEvent();
        end
        
        % Modify program settings based on key input
        if (~isempty(key))
            switch(key.charCode)
                case keyCodes.switchStepSize    % Switch step size mode
                    stepModePos = stepModePos + 1;
                    if stepModePos > length(stepModes)
                        stepModePos = 1;
                    end
                    % Number of beeps indicates new step size position
                    if ~silent
                        for i = 1:stepModePos
                            Snd('Play',sin(0:5000));
                        end
                    end
                    fprintf('User switched step size to %g\n',...
                        (stepModes(stepModePos)/ (adjustment_length - 1)));
                    
                case keyCodes.foundMatch  % Subject found a match
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    fprintf('User found match at %g test, %g primary\n',...
                        testScales(testPos), p1Scales(primaryPos));
                    % Save match
                    matches = [matches;...
                        [testScales(testPos), p1Scales(primaryPos)]];
                    matchPositions = [matchPositions; [testPos, primaryPos]];
                    save(fileLoc, 'matches', 'matchPositions', 'p1', 'p2',...
                        'test', 'cal','primarySpdsNominal',...
                        'primarySpdsPredicted', 'testSpdsNominal',...
                        'testSpdsPredicted', 'primaryStartStops',...
                        'testStartStops','whitePrimaries', 'whiteSettings',...
                        'whiteStarts', 'whiteStops','whiteSpdNominal',...
                        'subjectID', 'sessionNum','annulusData','sInterval',...
                        'lInterval', 'adjustment_length', 'foveal',...
                        'observer', 'subjectSettings');
                    % Set step size to largest and reset lights (random
                    % position for live experiment, initial values for
                    % simulation)
                    stepModePos = 1;
                    if simulate
                        primaryPos = 1;
                        testPos = 1;
                    else
                        primaryPos = randi(adjustment_length);
                        testPos = randi(adjustment_length);
                    end
                    
                case keyCodes.quit % Quit
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    fprintf('User exited program\n');
                    stillLooping = false;
                    
                case keyCodes.switchLightOrder % Switch primary/test order
                    rev = ~rev;
                    if ~silent
                        % One beep = primary first, two beeps = test first
                        Snd('Play',sin(0:5000)/ 20);
                        if rev
                            Snd('Play',sin(0:5000)/ 20);
                        end
                    end
                    fprintf('User switched order of primary and test lights\n');
                    
                case keyCodes.idealMatch % Switch to showing ideal match
                    ideal = ~ideal;
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    if ideal
                        fprintf('User switched to ideal match\n');
                    else
                        fprintf('User switched off ideal match\n');
                    end
                    
                case keyCodes.increaseIntensity % Scale up test intensity
                    testPos = testPos + stepModes(stepModePos);
                    if testPos > adjustment_length
                        testPos = adjustment_length;
                        if ~silent
                            Snd('Play',sin(0:5000));
                        end
                        fprintf('User reached upper test limit\n');
                    end
                    if ~silent
                        Snd('Play',sin(0:5000)/100);
                    end
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case keyCodes.decreaseIntensity % Scale down test intensity
                    testPos = testPos - stepModes(stepModePos);
                    if testPos < 1
                        testPos = 1;
                        if ~silent
                            Snd('Play',sin(0:5000));
                        end
                        fprintf('User reached lower test limit\n');
                    end
                    if ~silent
                        Snd('Play',sin(0:5000)/100);
                    end
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case keyCodes.increaseP1 % Move towards p1
                    primaryPos = primaryPos + stepModes(stepModePos);
                    if primaryPos > adjustment_length
                        primaryPos = adjustment_length;
                        if ~silent
                            Snd('Play',sin(0:5000));
                        end
                        fprintf('User reached upper primary limit\n');
                    end
                    if ~silent
                        Snd('Play',sin(0:5000)/100);
                    end
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    
                case keyCodes.decreaseP1 % Move towards p2
                    primaryPos = primaryPos - stepModes(stepModePos);
                    if primaryPos < 1
                        primaryPos = 1;
                        if ~silent
                            Snd('Play',sin(0:5000));
                        end
                        fprintf('User reached lower primary limit\n');
                    end
                    if ~silent
                        Snd('Play',sin(0:5000)/100);
                    end
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                        testScales(testPos), p1Scales(primaryPos));
            end
        end
    end
    % Once the time has elapsed, switch to the next light to display
    lightModePos = lightModePos + 1;
    if lightModePos > length(lights)
        lightModePos = 1;
    end
end

%% Close up
ol.setAll(false);
if simulate
    ListenChar(0);
else
    if ~foveal
        GLW_CloseAnnularStimulus();
    end
end

%% Optional Plots
% Optional plots
plot_responses = true;
if plot_responses
    % Separate subplots of primary/test trajectories
    figure();
    hold on; 
    subplot(2,1,1);
    plot(subjectSettings(:,2),'b*-');
    yline(subjectSettings(end,2)); 
    plot(length(subjectSettings(:,2)) + 1, p1Scales(pIdealIndex), 'r*');
    title('Primary Ratio');
    xlabel('Trial');
    ylabel('Proportion Red (p1)');
    
    subplot(2,1,2);
    plot(subjectSettings(:,1), 'b*-');
    yline(subjectSettings(end,1));
    plot(length(subjectSettings(:,1)) + 1, testScales(tIdealIndex), 'r*');  
    title('Test Intensity');
    xlabel('Trial');
    ylabel('Proportion of Maximal Intensity');
    sgtitle('Subject Settings Over Time');
    
    % Plot with both trajectories in tandem
    figure();
    hold on; 
    plot(subjectSettings(:,1), subjectSettings(:,2), 'b. ');
    plot(testScales(tIdealIndex), p1Scales(pIdealIndex), 'r*'); 
    xlabel('Proportional Test Intensity');
    ylabel('Primary Ratio (Proportion Red)');
    title('Subject Settings Over Time');
    legend('Subject Settings', 'Nominal Match'); 
end
end