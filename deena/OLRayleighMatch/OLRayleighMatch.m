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
%    'p1'             - integer wavelength of the first primary light in
%                       nm. Default is 670.
%    'p2'             - integer wavelength of the second primary light in
%                       nm. Default is 560.
%    'test'           - integer wavelength of the test light in nm. Default
%                       is 580.
%    'p1Scale'        - Numerical scale factor for the first primary light,
%                       between 0 and 1. Default is 1.
%    'p2Scale'        - Numerical scale factor for the second primary
%                       light, between 0 and 1. Default is 0.02.
%    'testScale'      - Numerical scale factor for the test light, between
%                       0 and 1. Default is 0.07.
%    'sInterval'      - length of time in s that the short stimulus is
%                       displayed for (or both stimuli when white light is
%                       used). Default is 0.25.
%    'lInterval'      - length of time in s that the long stimulus is
%                       displayed for (or the iti white light if it is
%                       used). Default is 1.
%    'foveal'         - logical indicating whether we are making foveal
%                       matches, in which case the annulus is not turned
%                       on. Default is true.
%    'white'          - logical indicating to run  a version of the
%                       experiment where the lights  are displayed with
%                       white in between. Default is false.
%    'simulate'       - logical. Set to true to run in simulated mode.
%                       Default is true.
%    'silent'         - logical. Set to true to turn off feedback beeps.
%                       Default is false.
%    'switchInterval' - When using a simulated observer, number of
%                       adjustments to make in a particular dimension
%                       (luminance or RG) before switching. Default is 1.
%    'numReversals'   - When using a simulated observer, number of
%                       reversals required before changing step size. Enter
%                       as a 2-element vector - the first element is the
%                       number of reversals for intermediate step sizes,
%                       the second is the number needed for the smallest
%                       step size. Default is [1 4].
%    'nBelowThreshold'      - When using a simulated observer with 
%                             threshold matching, number of pairs below 
%                             threshold required before recording a match. 
%                             Default is 1.
%    'thresholdScaleFactor' - When using a simulated observer with
%                             threshold matching, scale factor for matching
%                             threshold. Default is 2. 

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
%   06/04/20  dce       Added plotting option
%   06/09/20  dce       Added primary/test scale factors and threshold
%                       matching

%% Close any stray figures
close all;

%% Parse input
% This currently sets the primary and test wavelengths.
p = inputParser;
p.addParameter('p1', 670, @(x) (isnumeric(x)));
p.addParameter('p2', 560, @(x) (isnumeric(x)));
p.addParameter('test', 600, @(x) (isnumeric(x)));
p.addParameter('p1Scale', 1, @(x) (isnumeric(x)));
p.addParameter('p2Scale', 0.02, @(x) (isnumeric(x)));
p.addParameter('testScale', 0.07, @(x) (isnumeric(x)));
p.addParameter('sInterval', 0.25, @(x) (isnumeric(x)));
p.addParameter('lInterval', 1, @(x) (isnumeric(x)));
p.addParameter('foveal', true, @(x) (islogical(x)));
p.addParameter('white', false, @(x) (islogical(x)));
p.addParameter('simulate', true, @(x) (islogical(x)));
p.addParameter('silent', false, @(x) (islogical(x)));
p.addParameter('switchInterval', 1, @(x) (isnumeric(x)));
p.addParameter('numReversals', [1 4], @(x) (isnumeric(x)));
p.addParameter('nBelowThreshold', 1, @(x) (isnumeric(x)));
p.addParameter('thresholdScaleFactor', 2, @(x) (isnumeric(x)));

p.parse(varargin{:});
p1 = p.Results.p1;
p2 = p.Results.p2;
p1Scale = p.Results.p1Scale;
p2Scale = p.Results.p2Scale;
testScale = p.Results.testScale;
test = p.Results.test;
sInterval = p.Results.sInterval;
lInterval = p.Results.lInterval;
foveal = p.Results.foveal;
white = p.Results.white;
simulate = p.Results.simulate;
silent = p.Results.silent;
switchInterval = p.Results.switchInterval;
numReversals = p.Results.numReversals;
if (length(numReversals) ~= 2)
    error('Reversal vector must be 2x1');
end
nBelowThreshold = p.Results.nBelowThreshold;
thresholdScaleFactor = p.Results.thresholdScaleFactor; 

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
observer = [];
noisy = false;
thresholdMatching = false;
if simulate
    fprintf('Simulation Type:\n')
    fprintf('[1]: Simulated observer - threshold matches ("adjustment")\n');
    fprintf('[2]: Simulated observer - reversal matches ("forced choice")\n');
    fprintf('[3]: Do not use simulated observer\n');
    res = GetInput('Select option', 'number', 1);
    if res == 1 || res == 2
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
        fprintf('Add noise to decisions?\n')
        fprintf('[1]: Yes\n');
        fprintf('[2]: No\n');
        n = GetInput('Select option', 'number', 1);
        if n == 1
            noisy = true;
        end
        if res == 1
            thresholdMatching = true;   % We are making threshold matches
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
file = sprintf('OLRayleighMatchFineSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
    p1, p2, test, p1Scale, p2Scale, testScale);
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', file);
if ~exist(fName, 'file')
    OLRayleighMatchLightSettings(p1, p2, test, 'p1ScaleFactor', p1Scale,...
        'p2ScaleFactor', p2Scale, 'testScaleFactor', testScale);
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
adjustment_length = lightSettings.adjustmentLength;

%% Intialize OneLight and button box
ol = OneLight('simulate', simulate, 'plotWhenSimulating', false);
if simulate && ~sim_observer
    ListenChar(2);
    FlushEvents;
elseif ~simulate
    gamePad = GamePad();
end

%% Set up projector (if not making foveal matches)
annulusData = []; % Empty placeholder variable
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
            GLW_AnnularStimulusButtonBox();
        end
        % Run program to reset the annulus if its file does not exist
    else
        GLW_AnnularStimulusButtonBox();
    end
    % Load and display annulus file
    annulusData = load(annulusFile);
    annulusData.win.open;
    annulusData.win.draw;
    fprintf('\nProjector ready. Starting display loop\n')
end

%% Optional Plots
% If plot_responses is true, a pair of plots will be made for each match
plot_responses = true;
if plot_responses
    % Arrays for storing plot data
    nPlots = 2;           % Counter for figure index of current match
    matchSettingInd = 1;  % Position in subjectSettings where the current match starts
    
    % Separate subplots of primary/test trajectories
    figure(2);
    sgtitle('Subject Settings Over Time, Match 1');
    subplot(2,1,1);
    hold on;
    title('Primary Ratio');
    ylim([0 1.3]);
    xlabel('Number Adjustment');
    ylabel('Proportion Red (p1)');
    
    subplot(2,1,2);
    hold on;
    title('Test Intensity');
    ylim([0 1.3]);
    xlabel('Number Adjustment');
    ylabel('Proportion of Maximal Intensity');
    
    % Plot with both trajectories in tandem.
    figure(3);
    hold on;
    xlim([0 1]);
    ylim([0 1]);
    xlabel('Proportional Test Intensity');
    ylabel('Primary Ratio (Proportion Red)');
    title('Subject Settings Over Time, Match 1');
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
% script.
if foveal
    pIdealIndex = 116;
    tIdealIndex = 91;
else
    pIdealIndex = 125;
    tIdealIndex = 89;
end

%% Setup for display loop
% Initial settings
primaryPos = 1;      % Start with first position in primary array
testPos = 1;         % Start with first position in test array
stepModePos = 1;     % Start with the largest step size
lightModePos = 1;    % Start with the first light in the lights array
rev = false;         % Start with forward, not reverse order
ideal = false;       % Do not start with the ideal match
stillLooping = true; % Start looping

% Data-storing arrays
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array
subjectSettings = [testScales(testPos), p1Scales(primaryPos)];

% Settings for simulated observer
if sim_observer
    adjustingP = true;    % Start by adjusting the primary ratio
    adjustmentCount = 1;  % How long you've been adjusting one parameter
    
    p1_up_prev = true;    % Indicate that the program directed p1 and t to
    t_up_prev = true;     % increase on prior trials.
    
    pStepPos = 1;         % Start both primary and test adjustment at
    tStepPos = 1;         % largest step size
    
    numReversalsP = 0;    % Initially, have made 0 reversals for each
    numReversalsT = 0;    % adjustment
    if thresholdMatching
        % Number of consecutive light pairs below the matching threshold
        countBelowThreshold = 0;
    end
end

%% Display loop
% Loop through primary and test light until the user presses a key or the
% simulated observer decides to simulate a keypress
while(stillLooping)
    % Determine which light indices we are using
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
        % Simulated observer - perform calculations to figure out which
        % adjustment to make
        if sim_observer
            % Prompt the observer for a decision
            [p1_up, t_up, isBelowThreshold] = ...
                observerRayleighDecision(observer,primarySpdsPredicted(:, primaryPos),...
                testSpdsPredicted(:, testPos), 'noisy', noisy,...
                'thresholdScale', thresholdScaleFactor);
            % If we are making threshold matches, record whether this
            % combination of lights is below the threshold
            if thresholdMatching
                if isBelowThreshold
                    countBelowThreshold = countBelowThreshold + 1;
                    fprintf('Below threshold: %g\n', countBelowThreshold); 
                else
                    countBelowThreshold = 0;
                end
            end
            % Quit if all matches have been made
            [row, ~] = size(matches);
            if row == nMatches
                key.charCode = 'q';
            elseif (~thresholdMatching && pStepPos == length(stepModes)...
                    && tStepPos == length(stepModes)...
                    && (numReversalsP >= numReversals(2))...
                    && (numReversalsT >= numReversals(2)))
                % In the forced-choice version of the simulation, record
                % match and reset loop if both adjustments have had the
                % required number of reversals at the smallest step size
                key.charCode = ' ';
                p1_up_prev = logical(round(rand()));
                t_up_prev = logical(round(rand()));
                adjustingP = logical(round(rand())); 
                adjustmentCount = 1;
                pStepPos = 1;
                tStepPos = 1;
                numReversalsP = 0;
                numReversalsT = 0;
            elseif thresholdMatching &&...
                    countBelowThreshold == nBelowThreshold
                % In the adjustment version of the simulation, record match
                % and reset loop if the required number of matches has been
                % made
                key.charCode = ' ';
                p1_up_prev = logical(round(rand()));
                t_up_prev = logical(round(rand()));
                adjustingP = logical(round(rand()));
                adjustmentCount = 1;
                pStepPos = 1;
                tStepPos = 1;
                numReversalsP = 0;
                numReversalsT = 0;
                countBelowThreshold = 0;
            elseif adjustingP  % Adjust primary if match is not yet reached
                stepModePos = pStepPos;
                if p1_up == p1_up_prev % Continue in existing direction
                    if p1_up
                        key.charCode = 'r';
                    else
                        key.charCode = 'l';
                    end
                    adjustmentCount = adjustmentCount + 1;
                else        % Reversal
                    % Move as indicated, and lower step size if the number
                    % of reversals for that stage has been reached.
                    if numReversalsP < numReversals(1) || stepModePos == length(stepModes)
                        if p1_up
                            key.charCode = 'r';
                        else
                            key.charCode = 'l';
                        end
                        adjustmentCount = adjustmentCount + 1;
                        numReversalsP = numReversalsP + 1;
                        fprintf('Primary reversal %g, step size %g\n',...
                            numReversalsP, (stepModes(stepModePos)/ (adjustment_length - 1)));
                    else
                        key.charCode = 's';
                        pStepPos = pStepPos + 1;
                        numReversalsP = 0;
                    end
                end
                p1_up_prev = p1_up;
            else   % Adjust test intensity if match is not yet reached
                stepModePos = tStepPos;
                if t_up == t_up_prev % Continue in existing direction
                    if t_up
                        key.charCode = 'u';
                    else
                        key.charCode = 'd';
                    end
                    adjustmentCount = adjustmentCount + 1;
                else        % Reversal
                    % Move as indicated, and lower step size if the number
                    % of reversals for that stage has been reached.
                    if stepModePos == length(stepModes) || numReversalsT < numReversals(1)
                        if t_up
                            key.charCode = 'u';
                        else
                            key.charCode = 'd';
                        end
                        numReversalsT = numReversalsT + 1;
                        fprintf('Test reversal %g, step size %g\n',...
                            numReversalsT, (stepModes(stepModePos)/ (adjustment_length - 1)));
                        adjustmentCount = adjustmentCount + 1;
                    else
                        key.charCode = 's'; % Lower step size
                        tStepPos = tStepPos + 1;
                        numReversalsT = 0;
                    end
                end
                t_up_prev = t_up;
            end
            % Switch which light you're adjusting if the limit is reached
            if adjustmentCount > switchInterval
                adjustingP = ~adjustingP;
                adjustmentCount = 1;
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
                    if sim_observer && adjustingP
                        fprintf('User switched primary step size to %g\n',...
                            (stepModes(stepModePos)/ (adjustment_length - 1)));
                    elseif sim_observer
                        fprintf('User switched test step size to %g\n',...
                            (stepModes(stepModePos)/ (adjustment_length - 1)));
                    else
                        fprintf('User switched step size to %g\n',...
                            (stepModes(stepModePos)/ (adjustment_length - 1)));
                    end
                    
                case keyCodes.foundMatch  % Subject found a match
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    % Calculate match. If this is a live experiment or a
                    % simulated threshold match, simply provide the subject
                    % setting. For a simulated observer with a
                    % forced-choice match, average the last two distinct
                    % light settings.
                    if sim_observer && ~thresholdMatching
                        % Find previous values to average with
                        pUnique = unique(flip(subjectSettings(:, 2)), 'stable');
                        tUnique = unique(flip(subjectSettings(:, 1)), 'stable');
                        pMatch = mean([pUnique(1), pUnique(2)]);
                        tMatch = mean([tUnique(1), tUnique(2)]);
                    else
                        pMatch = p1Scales(primaryPos);
                        tMatch = testScales(testPos);
                    end
                    fprintf('User found match at %g test, %g primary\n',...
                        tMatch, pMatch);
                    matches = [matches; [tMatch, pMatch]];
                    matchPositions = [matchPositions; [testPos, primaryPos]];
                    save(fileLoc, 'matches', 'matchPositions',...
                        'subjectSettings', 'p1', 'p2', 'test', 'cal',...
                        'primarySpdsNominal','primarySpdsPredicted',...
                        'testSpdsNominal', 'testSpdsPredicted',...
                        'primaryStartStops','testStartStops','subjectID',...
                        'sessionNum','annulusData','sInterval','lInterval',...
                        'adjustment_length', 'foveal', 'white', 'simulate',...
                        'numReversals', 'switchInterval', 'observer',...
                        'noisy', 'p1Scale', 'p2Scale', 'testScale',...
                        'thresholdMatching', 'thresholdScaleFactor',...
                        'nBelowThreshold');
                    if plot_responses
                        nAdjustments = length(subjectSettings(matchSettingInd:end,1));
                        matchSettingInd = length(subjectSettings(:,1))+1;
                        % Current individual trajectory figure
                        figure(nPlots);
                        subplot(2,1,1)
                        p1 = plot(nAdjustments+1, matches(end,2), 'r* ',...
                            nAdjustments+1, p1Scales(pIdealIndex), 'gs',...
                            'MarkerSize', 7);
                        legend(p1,'Subject Match','Nominal Match');
                        subplot(2,1,2);
                        p2 = plot(nAdjustments+1, matches(end,1), 'r* ',...
                            nAdjustments+1, testScales(tIdealIndex), 'gs',...
                            'MarkerSize', 7);
                        legend(p2,'Subject Match','Nominal Match');
                        
                        % Current joint trajectory figure
                        figure(nPlots+1);
                        p3 = plot(matches(end,1), matches(end,2), 'r* ', 'MarkerSize', 10, 'LineWidth', 1.5);
                        p4 = plot(testScales(tIdealIndex), p1Scales(pIdealIndex), 'gs',...
                            'MarkerFaceColor', 'g');
                        legend([p3,p4], 'Subject Match', 'Nominal Match');
                        
                        % Add new figures for the next match
                        nPlots = nPlots + 2;
                        figure(nPlots);
                        title1 = sprintf('Subject Settings Over Time, Match %g', nPlots/2);
                        sgtitle(title1);
                        subplot(2,1,1);
                        hold on;
                        title('Primary Ratio');
                        ylim([0 1.3]);
                        xlabel('Number Adjustment');
                        ylabel('Proportion Red (p1)');
                        
                        subplot(2,1,2);
                        hold on;
                        title('Test Intensity');
                        ylim([0 1.3]);
                        xlabel('Number Adjustment');
                        ylabel('Proportion of Maximal Intensity');
                        
                        figure(nPlots +1);
                        hold on;
                        xlim([0 1]);
                        ylim([0 1]);
                        xlabel('Proportional Test Intensity');
                        ylabel('Primary Ratio (Proportion Red)');
                        title2 = sprintf('Subject Settings Over Time, Match %g', nPlots/2);
                        title(title2);
                    end
                    % Set step size to largest and reset lights (random
                    % position for live experiment, initial values for
                    % simulation)
                    stepModePos = 1;
                    primaryPos = randi(adjustment_length);
                    testPos = randi(adjustment_length);
                    
                case keyCodes.quit % Quit
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    fprintf('User exited program\n');
                    stillLooping = false;
                    break;
                    
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
                    fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                        testScales(testPos), p1Scales(primaryPos));
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    if plot_responses
                        figure(nPlots);
                        subplot(2,1,1);
                        plot(subjectSettings(matchSettingInd:end,2),'b*-');
                        subplot(2,1,2);
                        plot(subjectSettings(matchSettingInd:end,1),'b*-');
                        figure(nPlots+1);
                        plot(subjectSettings(matchSettingInd:end,1),...
                            subjectSettings(matchSettingInd:end,2), 'b*-');
                    end
                    
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
                    fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                        testScales(testPos), p1Scales(primaryPos));
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    if plot_responses
                        figure(nPlots);
                        subplot(2,1,1);
                        plot(subjectSettings(matchSettingInd:end,2),'b*-');
                        subplot(2,1,2);
                        plot(subjectSettings(matchSettingInd:end,1),'b*-');
                        figure(nPlots+1);
                        plot(subjectSettings(matchSettingInd:end,1),...
                            subjectSettings(matchSettingInd:end,2), 'b*-');
                    end
                    
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
                    fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                        testScales(testPos), p1Scales(primaryPos));
                    subjectSettings = [subjectSettings; [testScales(testPos), p1Scales(primaryPos)]];
                    if plot_responses
                        figure(nPlots);
                        subplot(2,1,1);
                        plot(subjectSettings(matchSettingInd:end,2),'b*-');
                        subplot(2,1,2);
                        plot(subjectSettings(matchSettingInd:end,1),'b*-');
                        figure(nPlots+1);
                        plot(subjectSettings(matchSettingInd:end,1),...
                            subjectSettings(matchSettingInd:end,2), 'b*-');
                    end
                    
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
                    if plot_responses
                        figure(nPlots);
                        subplot(2,1,1);
                        plot(subjectSettings(matchSettingInd:end,2),'b*-');
                        subplot(2,1,2);
                        plot(subjectSettings(matchSettingInd:end,1),'b*-');
                        figure(nPlots+1);
                        plot(subjectSettings(matchSettingInd:end,1),...
                            subjectSettings(matchSettingInd:end,2), 'b*-');
                    end
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
elseif ~foveal
    GLW_CloseAnnularStimulus();
end

% Close extra plots
if plot_responses && sim_observer
    close(figure(nPlots), figure(nPlots+1));
end
end