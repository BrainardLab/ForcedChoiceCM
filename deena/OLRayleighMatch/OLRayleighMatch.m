function OLRayleighMatch(subjectID,sessionNum,varargin)
% Program to run a Rayleigh match experiment on the OneLight, or as a
% simulation
%
% Syntax:
%   OLRayleighMatch(subjectID, sessionNum)
%
% Description:
%    Displays a mixture of two primary lights followed by a test light in
%    alternation.
%
%    Default wavelengths are 600 for test, 670 for p1, and 560 for
%    p2, but these can also be entered by the user as key/value pairs.
%
%    The live experiment can be ran in two ways: "adjustment" or "forced
%    choice." Under the forced choice approach, subjects are prompted to
%    adjust the primary ratio and the test light intensity using the game
%    pad's directional pad (right moves the primary ratio towards p1 and
%    left towards p2, up increases the test intensity and down decreases
%    the test intensity.) However, the program determines when to change
%    the step size and when to declare a match. The forced-choice matches
%    are an average of the last four settings.
%
%    In the "adjustment" procedure, the subject can use the top 'Y' button
%    to toggle primary step size, the left 'X' button to toggle test step
%    size, the right 'B' button to record a match, the 'Back' button to
%    show the ideal match, and the 'Start' button for switching light
%    order, in addition to the controls described above. In both
%    procedures, the bottom 'A' button can be used to force quit.
%
%    In addition to the live experiment, the program can run a simulation
%    by creating a simulated observer based on the Asano model and searing
%    for the best match for this observer. The experimenter can choose
%    between two decision rules - a "forced choice" procedure that
%    records a match based on the number of reversals, and an "adjustment"
%    procedure that records a match if the two lights meet a similarity
%    threshold.
%
% Inputs:
%    subjectID        - Character vector of subject ID
%    sessionNum       - Integer session number.
%
% Outputs:
%    Saves a .mat file titled by subject and session number, which includes
%    a record of subjects' matches and various experimental parameters.
%    There is also a setting to save figures as well.
%
% Optional key-value pairs:
%    'p1'             - integer wavelength of the first primary light in
%                       nm. Default is 670.
%    'p2'             - integer wavelength of the second primary light in
%                       nm. Default is 560.
%    'test'           - integer wavelength of the test light in nm. Default
%                       is 600.
%    'p1Scale'        - Numerical scale factor for the first primary light,
%                       between 0 and 1. Default is 1.
%    'p2Scale'        - Numerical scale factor for the second primary
%                       light, between 0 and 1. Default is 0.02.
%    'testScale'      - Numerical scale factor for the test light, between
%                       0 and 1. Default is 0.1.
%    'sInterval'      - length of time in s that each light is
%                       displayed for. Default is 0.5.
%    'lInterval'      - length of time in s between lights. Default is 1.
%    'fieldSize'         - logical indicating whether we are making foveal
%                       matches, in which case the annulus is not turned
%                       on. Default is true.
%    'age'            - Integer age for observer. Default is 32.
%    'resetAnnulus'   - logical indicating to run a script that lets the
%                       experimenter reset the annulus. Default is false.
%    'plotResponses'  - Logical indicating to make plots of the
%                       experimental timecourse and print output. Default
%                       is true.
%    'savePlots'      - Logical indicating to save figures. Default is
%                       true.
%    'silent'         - logical. Set to true to turn off feedback beeps.
%                       Default is true.
%    'simObserver'    - logical. Set to true to run with a simulated
%                       observer. Default is true.
%    'observerParams' - 1x8 vector of Asano individual difference params
%                       for simulated observer. Default is zeros(1,8)
%    'opponentParams' - 1x4 vector of opponent contrast parameters. Default
%                       is [40.3908  205.7353   62.9590    1.0000].
%    'nReversals'     - When using a simulated observer, number of
%                       reversals required before changing step size. Enter
%                       as a 2-element vector - the first element is the
%                       number of reversals for intermediate step sizes,
%                       the second is the number needed for the smallest
%                       step size. Default is [1 4].
%    'adjustment'     - Logical indicating to use an adjustment method, not
%                       forced choice. For simulated observers, this means
%                       that matches are recorded when the cone excitation
%                       difference for the 2 spectra falls below a
%                       threshold. For live observers, this means that
%                       observers can freely record step sizes and make a
%                       match. Default is false.
%    'nBelowThreshold'      - When using a simulated observer with
%                             adjustment matching, number of pairs below
%                             threshold required before recording a match.
%                             Default is 1.
%    'thresholdScaleFactor' - When using a simulated observer with
%                             adjustment matching, scale factor for matching
%                             threshold. Default is 0.5.
%     'noiseScaleFactor'     -Number >=0 which determines which scalar
%                             multiple of the opponent noise SD should be
%                             used as the observer noise SD. Default is 0.
%    'nObserverMatches'     - Number of matches to simulate/run. Default is
%                             1.
%    'adjustmentLength'     - Number of possible steps available for
%                             adjusting the primary and test lights.
%                             Default is 201.
%    'simNominalLights'     - Logical indicating to run the simulation with
%                             nominal, rather than predicted, spds. Default
%                             is false.
%    'monochromatic'        - Logical indicating to run the simulation with
%                             monochromatic spds. Default is false.
%    'monochromaticS'       - Vector with wavelength sampling if
%                             monochromatic lights are used. Default is
%                             [380 2 201].
%    'lambdaRef'           - Number between 0 and 1 indicating which value
%                            of lambda to use for a reference primary when
%                            calculating simulated opponent contrasts. Must
%                            be a member of p1Scales. When empty, opponent
%                            contrasts are not computed relative to a
%                            reference. Default is [].
%    'stimLimits'           -4-element vector to store limits on stimulus
%                            parameters. Follows the format [min lambda,
%                            max lambda, min test intensity, max test
%                            intensity].Each element must be between 0 and
%                            1. Default is [].
%    'refFirst'             -Logical. If true, displays the test light
%                            before the primary mixture (and for a longer
%                            period of time, instead of the other way
%                            around. Default is false. 
%    'pairStepSizes'        -Logical. If true, adjusts primary and test
%                            step sizes together instead of separately. 
%                            Default is false.

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
%   06/15/20  dce       Removed prompts for user input
%   06/18/20  dce       Debugged plotting and nominal match setting
%   07/06/20  dce       Switched nominal match calculation to analytical
%                       method.
%   07/14/20  dce       Changed nominal match to use predicted spds.
%                       Changed default adjustment length.
%   07/21/20  dce       Changed method for calculating nominal matches
%   07/24/20  dce       Added exit condition for infinite looping in
%                       threshold case
%   08/05/20  dce       Added opponent contrast params input
%   08/09/20  dce       Changed ideal match to best available, not nominal
%   10/28/20  dce       Added monochromatic simulation matching
%   11/15/20  dce       Added optional stimulus limits, reference light
%   02/05/20  dce       Fixed charcode system for using real OneLight
%   03/07/20  dce       Restructured to incorporate a forced-choice
%                       procedure and make judgements for both lights on
%                       each trial
%   03/08/20  dce       Added flicker back in while waiting for subject
%                       decisions
%   03/15/20  dce       Edited for style, added two new key-value pairs


%% Close any stray figures
close all;

%% Parse input
p = inputParser;
p.addParameter('p1',670,@(x)(isnumeric(x)));
p.addParameter('p2',560,@(x)(isnumeric(x)));
p.addParameter('test',600,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.1,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',201,@(x)(isnumeric(x)));
p.addParameter('sInterval',0.5,@(x)(isnumeric(x)));
p.addParameter('lInterval',1,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('resetAnnulus',false,@(x)(islogical(x)));
p.addParameter('plotResponses',true,@(x) (islogical(x)));
p.addParameter('savePlots',true,@(x) (islogical(x)));
p.addParameter('silent',true,@(x)(islogical(x)));
p.addParameter('simObserver',true,@(x)(islogical(x)));
p.addParameter('observerParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('adjustment',false,@(x)(islogical(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('noiseScaleFactor',0,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('simNominalLights',false,@(x)(islogical(x)));
p.addParameter('monochromatic',false,@(x)(islogical(x)));
p.addParameter('monochromaticS',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('lambdaRef',[],@(x)(isnumeric(x)));
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.addParameter('testFirst',false,@(x)(islogical(x)));
p.addParameter('pairStepSizes',false,@(x)(islogical(x)));
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
p1Scale = p.Results.p1Scale;
p2Scale = p.Results.p2Scale;
testScale = p.Results.testScale;
adjustmentLength = p.Results.adjustmentLength;
sInterval = p.Results.sInterval;
lInterval = p.Results.lInterval;
fieldSize = p.Results.fieldSize;
age = p.Results.age;
resetAnnulus = p.Results.resetAnnulus;
silent = p.Results.silent;
simObserver = p.Results.simObserver;
observerParams = p.Results.observerParams;
opponentParams = p.Results.opponentParams;
nObserverMatches = p.Results.nObserverMatches;
nReversals = p.Results.nReversals;
adjustment = p.Results.adjustment;
nBelowThreshold = p.Results.nBelowThreshold;
thresholdScaleFactor = p.Results.thresholdScaleFactor;
noiseScaleFactor = p.Results.noiseScaleFactor;
plotResponses = p.Results.plotResponses;
savePlots = p.Results.savePlots;
simNominalLights = p.Results.simNominalLights;
monochromatic = p.Results.monochromatic;
monochromaticS = p.Results.monochromaticS;
stimLimits = p.Results.stimLimits;
lambdaRef = p.Results.lambdaRef;
testFirst = p.Results.testFirst;
pairStepSizes = p.Results.pairStepSizes;

% Input error checking
if (length(nReversals)~=2)
    error('Reversal vector must be 2x1');
end
if (length(observerParams)~=8)
    error('Observer parameters must be entered as an 8-element vector');
end
if (length(opponentParams)~=4)
    error('Opponent contrast parameters must be entered as an 4-element vector');
end
if (monochromatic && ~simObserver)
    error('Monochromatic matching can only be used with a simulation method');
end

%% Set up directory for saving results
% Create directory named subjectID for saving data, if it doesn't exist
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Create data file with name subject ID_SessionNumber. Throw error if a
% file already exists with that name.
fileName = [subjectID,'_',num2str(sessionNum),'.mat'];
fileLoc = fullfile(outputDir,fileName);
if exist(fileLoc,'file')
    error('Specified output file %s already exists',fileName);
end

%% Find light settings
% Find precomputed spectra, or compute if they do not exist
calName = [];
lightFileName = [];
if ~monochromatic  % Matching using OL primaries
    lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
        adjustmentLength,p1,p2,test,p1Scale,p2Scale,testScale);
    lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops',lightFile);
    if ~exist(lightFileName,'file')
        OLRayleighMatchLightSettings(p1,p2,test,'p1ScaleFactor',p1Scale,...
            'p2ScaleFactor',p2Scale,'testScaleFactor',testScale,...
            'adjustmentLength',adjustmentLength);
    end
    % Save some variables locally
    lightSettings = load(lightFileName);
    calName = lightSettings.cal.describe.calType;
    primarySpdsNominal = lightSettings.primarySpdsNominal;
    primarySpdsPredicted = lightSettings.primarySpdsPredicted;
    testSpdsNominal = lightSettings.testSpdsNominal;
    testSpdsPredicted = lightSettings.testSpdsPredicted;
    primaryStartStops = lightSettings.primaryStartStops;
    testStartStops = lightSettings.testStartStops;
    p1Scales = lightSettings.p1Scales;
    testScales = lightSettings.testScales;
    S = lightSettings.cal.computed.pr650S;
    
else  % Monochromatic matching
    p1Scales = linspace(0,1,adjustmentLength);
    p2Scales = 1-p1Scales;
    testScales = linspace(0,1,adjustmentLength);
    
    S = monochromaticS;
    wls = SToWls(S);
    primarySpdsNominal = zeros(length(wls),length(p1Scales));
    primarySpdsNominal(wls==p1,:) = 1*p1Scale*p1Scales;
    primarySpdsNominal(wls==p2,:) = 1*p2Scale*p2Scales;
    primarySpdsPredicted = primarySpdsNominal;
    
    testSpdsNominal = zeros(length(wls),length(testScales));
    testSpdsNominal(wls==test,:) = 1*testScale*testScales;
    testSpdsPredicted = testSpdsNominal;
end

%% Set up limits on stimulus parameters. By default, all available indices
% are allowed
allowedLambdaInds = 1:adjustmentLength;
allowedTIInds = 1:adjustmentLength;
if ~isempty(stimLimits)
    allowedLambdaInds = allowedLambdaInds(p1Scales >= stimLimits(1)...
        & p1Scales <= stimLimits(2));
    allowedTIInds = allowedTIInds(testScales >= stimLimits(3)...
        & testScales <= stimLimits(4));
    
    % Fill in so there are at least 2 values in each vector
    if isempty(allowedLambdaInds) || length(allowedLambdaInds)< 4
        [~,minInd] = min(abs(p1Scales - stimLimits(1)));
        [~,maxInd] = min(abs(p1Scales - stimLimits(2)));
        allowedLambdaInds = unique([minInd-1,minInd,maxInd,maxInd+1]);
    end
    if isempty(allowedTIInds) || length(allowedTIInds) <4
        [~,minInd] = min(abs(testScales - stimLimits(3)));
        [~,maxInd] = min(abs(testScales - stimLimits(4)));
        allowedTIInds = unique([minInd-1,minInd,maxInd,maxInd+1]);
    end
end

%% Initialize simulated observer if desired
stdObserver = genRayleighObserver('fieldSize',fieldSize,'age',age,...
    'coneVec',zeros(1,8),'opponentParams',opponentParams,...
    'S',S);
if simObserver
    observer = genRayleighObserver('fieldSize',fieldSize,'age',age,...
        'coneVec',observerParams,'opponentParams',opponentParams,...
        'S',S);
else
    observer = stdObserver;
end

%% Find reference spd, if desired
refSpd = [];
% Which spds are we using in simulation?
if simNominalLights
    primarySpds = primarySpdsNominal;
    testSpds = testSpdsNominal;
else
    primarySpds = primarySpdsPredicted;
    testSpds = testSpdsPredicted;
end
if ~isempty(lambdaRef)
    refSpd = primarySpds(:,p1Scales==lambdaRef);
    if isempty(refSpd)
        error('Provided reference lambda is not a selected primary mixture scalar');
    end
end

%% Find nominal match
[idealTestSpd,idealPrimarySpd,tIdealIndex,pIdealIndex] =...
    searchPredictedRayleighMatch(testSpds,primarySpds,observer,'refSpd',refSpd);
idealTestIntensity = testScales(tIdealIndex);
idealPRatio = p1Scales(pIdealIndex);

%% Intialize OneLight and button box
if ~simObserver
    ol = OneLight();
    gamePad = GamePad();
end

%% Set up projector (if not making foveal or simulated matches)
annulusData = []; % Empty placeholder variable
% Check if an annulus file exists. If it does, can choose whether to use
% the existing file or reset the annulus.
if ~simObserver && fieldSize > 2
    annulusFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'projectorSettings','OLAnnulusSettings.mat');
    % Reset annulus if the file does not exist or if the experimenter chose
    % to
    if ~exist(annulusFile,'file') || resetAnnulus
        GLW_AnnularStimulusButtonBox();
    end
    % Load and display annulus file
    annulusData = load(annulusFile);
    annulusData.win.open;
    annulusData.win.draw;
end

%% Display parameters
% Possible step sizes (relative to adjustmentLength)
stepModes = [floor(adjustmentLength/5),floor(adjustmentLength/20),...
    floor(adjustmentLength/100),floor(adjustmentLength/200),...
    floor(adjustmentLength/400),floor(adjustmentLength/800),...
    floor(adjustmentLength/1600),floor(adjustmentLength/3200)];
stepModes = stepModes(stepModes ~= 0);

% Initial setup 
ideal = false;              % Do not start with the ideal match
firstAdjustment = true;     % This is the first adjustment of the match
stillLooping = true;        % Start looping
displayLoopCounter = 0;     % Count number of iterations
p1_up_prev = true;
t_up_prev = true;

% Data-storing arrays
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array
subjectSettings = [];
primaryRes = [];            % Fill with 1's for every "up" trial, 0's for every "down" trial
testRes = [];               % Fill with 1's for every "up" trial, 0's for every "down" trial
pRevIndices = [];           % Indices where a reversal occurred
tRevIndices = [];           % Indices where a reversal occurred
pStepSwitchInds = [];       % Indices where step size was changed
tStepSwitchInds = [];       % Indices where step size was changed

% Optional plotting setup
if plotResponses
    nPlots = -1;           % Counter for current figure index
    matchSettingInd = 1;  % SubjectSettings position for current match start
end

%% Display loop
while(stillLooping)
    % Reset loop control parameters
    match = false;
    quit = false;
    switchTStepSize = false;
    switchPStepSize = false;
    p1_up = false;
    p1_down = false;
    t_up = false;
    t_down = false;
    displayLoopCounter = displayLoopCounter + 1;
    
     
    % Setup if it's a new match
    if firstAdjustment
        pStepModePos = 1;     % Start with the largest step sizes
        tStepModePos = 1;
        nReversalsP = 0;      % Initially, have made 0 reversals for each
        nReversalsT = 0;      % adjustment
        countBelowThreshold = 0;    % Count pairs below match threshold
        primaryPos = randsample(allowedLambdaInds,1);
        testPos = randsample(allowedTIInds,1);
        if monochromatic
            primaryPos = max(primaryPos,2);
            testPos = max(testPos,2);
        end
        subjectSettings = [subjectSettings;...
            [testScales(testPos),p1Scales(primaryPos)]];
        
        if plotResponses
            nPlots = nPlots+2;
            % Plot with primary and test settings separately
            plot1 = figure(nPlots);
            title1 = sprintf('Subject Settings Over Time, Match %g',ceil(nPlots/2));
            %             sgtitle(title1);
            subplot(2,1,1);
            hold on;
            title('Primary Ratio');
            ylim([0 1.3]);
            xlabel('Number Adjustment');
            ylabel('Proportion Red (p1)');
            subplot(2,1,2);
            hold on;
            title('Reference Intensity');
            ylim([0 1.3]);
            xlabel('Number Adjustment');
            ylabel('Proportion of Maximal Intensity');
            
            % Plot with both trajectories in tandem
            plot2 = figure(nPlots+1);
            hold on;
            xlim([0 1]);
            ylim([0 1]);
            ylabel('Proportional Reference Intensity');
            xlabel('Primary Ratio (Proportion Red)');
            title2 = sprintf('Subject Settings Over Time, Match %g',ceil(nPlots/2));
            title(title2);
            plot(idealPRatio,idealTestIntensity,'gs',...
                'MarkerFaceColor','g');
        end
    end
    
    % Define light indices
    if ideal
        pI = pIdealIndex;
        tI = tIdealIndex;
    else
        pI = primaryPos;
        tI = testPos;
    end
    
    if testFirst
        starts = {squeeze(testStartStops(tI,1,:))',squeeze(primaryStartStops(pI,1,:))'};
        stops = {squeeze(testStartStops(tI,2,:))',squeeze(primaryStartStops(pI,2,:))'};
    else
        starts = {squeeze(primaryStartStops(pI,1,:))', squeeze(testStartStops(tI,1,:))'};
        stops = {squeeze(primaryStartStops(pI,2,:))', squeeze(testStartStops(tI,2,:))'};
    end
    intervals = [lInterval sInterval];
    
    % In a forced-choice live experiment, we show the specified lights then
    % pause for a decision period. In a simulated experiment, we prompt the
    % simulated observer for a decision.
    if (~simObserver)
        % Get the observer decision. If using the forced-choice method, we
        % get one R/G and one brightness decision. If using the adjustment
        % method, we wait for the first keypress.
        waitingForResponse = true;
        innerLoopCounter = 1;
        if adjustment
            while(waitingForResponse)
                nowTime = mglGetSecs;
                innerLoopCounter = innerLoopCounter+1;
                while(mglGetSecs < nowTime+intervals(mod(innerLoopCounter,2)+1))  % Flicker while waiting for response
                    ol.setMirrors(starts{mod(innerLoopCounter,2)+1},stops{mod(innerLoopCounter,2)+1});
                    key = gamePad.getKeyEvent();
                    if (~isempty(key))
                        switch(key.charCode)
                            case 'GP:Y'
                                switchPStepSize = true;
                            case 'GP:X'
                                switchTStepSize = true;
                            case 'GP:B'
                                match = true;
                            case 'GP:A'
                                quit = true;
                            case 'GP:Back'
                                ideal = ~ideal;
                                if ~silent
                                    Snd('Play',sin(0:5000));
                                end
                                if plotResponses
                                    if ideal
                                        fprintf('User switched to ideal match\n');
                                    else
                                        fprintf('User switched off ideal match\n');
                                    end
                                end
                            case 'GP:Start'
                                testFirst = ~testFirst;
                                if ~silent
                                    % One beep = primary first, two beeps = test first
                                    Snd('Play',sin(0:5000)/ 20);
                                    if testFirst
                                        Snd('Play',sin(0:5000)/ 20);
                                    end
                                end
                                if plotResponses
                                    fprintf('User switched order of primary and test lights\n');
                                end
                            case 'GP:North'
                                t_up = true;
                            case 'GP:South'
                                t_down = true;
                            case 'GP:East'
                                p1_up = true;
                            case 'GP:West'
                                p1_down = true;
                        end
                        waitingForResponse = false;
                        break;
                    end
                end
            end
            
        else   % Forced choice case
            % Show the two lights, pausing on the long one 
            for i = 1:3
                nowTime = mglGetSecs;
                innerLoopCounter = innerLoopCounter+1;
                while(mglGetSecs < nowTime+intervals(mod(innerLoopCounter,2)+1))  % Flicker while waiting for response
                    ol.setMirrors(starts{mod(innerLoopCounter,2)+1},stops{mod(innerLoopCounter,2)+1});
                end 
            end
            innerLoopCounter = 1;
            
            % Prompt for R/G decision
            if testFirst
                Speak('Short Redness?');
            else
                Speak('Long Redness?');
            end
            while(waitingForResponse)
                nowTime = mglGetSecs;
                innerLoopCounter = innerLoopCounter+1;
                while(mglGetSecs < nowTime+intervals(mod(innerLoopCounter,2)+1))  % Flicker while waiting for response
                    ol.setMirrors(starts{mod(innerLoopCounter,2)+1},stops{mod(innerLoopCounter,2)+1});
                    key = gamePad.getKeyEvent();
                    if (~isempty(key))
                        switch(key.charCode)
                            case 'GP:East'
                                p1_up = true;
                                waitingForResponse = false;
                                break;
                            case 'GP:West'
                                p1_down = true;
                                waitingForResponse = false;
                                break;
                            case 'GP:A' % Force quit
                                quit = true;
                                waitingForResponse = false;
                                break;
                            otherwise
                                Speak('Invalid key');
                        end
                    end
                end
            end
            if ~silent
                Snd('Play',sin(0:5000)/100);
            end
            innerLoopCounter = 1;
            
            % Prompt for ref decision
            if testFirst
                Speak('Long Brightness?');
            else
                Speak('Short Brightness?');
            end
            waitingForResponse = true;
            while(waitingForResponse)
                nowTime = mglGetSecs;
                innerLoopCounter = innerLoopCounter+1;
                while(mglGetSecs < nowTime+intervals(mod(innerLoopCounter,2)+1))  % Flicker while waiting for response
                    ol.setMirrors(starts{mod(innerLoopCounter,2)+1},stops{mod(innerLoopCounter,2)+1});
                    key = gamePad.getKeyEvent();
                    if (~isempty(key))
                        switch(key.charCode)
                            case 'GP:North'
                                t_up = true;
                                waitingForResponse = false;
                                break;
                            case 'GP:South'
                                t_down = true;
                                waitingForResponse = false;
                                break;
                            case 'GP:A'  % Force quit
                                quit = true;
                                waitingForResponse = false;
                                break;
                            otherwise
                                Speak('Invalid key');
                        end
                    end
                end
            end
        end
        if ~silent
            Snd('Play',sin(0:5000)/100);
        end
    else  % Prompt simulated observer for decisions
        [p1_up,t_up,isBelowThreshold] = ...
            observerRayleighDecision(observer,primarySpds(:,primaryPos),...
            testSpds(:,testPos),'thresholdScale',thresholdScaleFactor,...
            'noiseScale',noiseScaleFactor,'refSpd',refSpd);
        if ~p1_up
            p1_down = true;
        end
        if ~t_up
            t_down = true;
        end
        
        %  A couple methods specific for simulation with adjustment method.
        if adjustment && isBelowThreshold
            countBelowThreshold = countBelowThreshold+1;
        end
    end
    
    % Check for a primary reversal
    if ((p1_up && ~p1_up_prev) || (p1_down && p1_up_prev)) && ~firstAdjustment % Reversal
        nReversalsP = nReversalsP+1;  % Count the reversal
        pRevIndices = [pRevIndices,displayLoopCounter];
        if plotResponses
            fprintf('Primary reversal %g, step size %g\n',...
                nReversalsP,(stepModes(pStepModePos)/...
                (adjustmentLength-1)));
        end
    end
    
    % Check for a ref reversal
    if ((t_up && ~t_up_prev) || (t_down && t_up_prev)) && ~firstAdjustment % Reversal
        nReversalsT = nReversalsT+1;  % Count the reversal
        tRevIndices = [tRevIndices,displayLoopCounter];
        if plotResponses
            fprintf('Test reversal %g, step size %g\n',...
                nReversalsT,(stepModes(tStepModePos)/...
                (adjustmentLength-1)));
        end
    end
    
    % Process the observer responses, if not using the live adjustment
    % experiment
    if ~(adjustment && ~simObserver)
        % Check if we're at a match
        if (simObserver && adjustment && countBelowThreshold ==...
                nBelowThreshold) || (~adjustment ...
                && pStepModePos == length(stepModes)...
                && tStepModePos == length(stepModes)...
                && (nReversalsP >= nReversals(2))...
                && (nReversalsT >= nReversals(2)))
            match = true;
        end
        
        % Check if primary step size needs to be adjusted
        if ((p1_up && ~p1_up_prev) || (p1_down && p1_up_prev))...
                && ~firstAdjustment && nReversalsP >= nReversals(1)...
                && pStepModePos ~= length(stepModes) 
                switchPStepSize = true;
        end
        
        % Check if test step size needs to be adjusted
        if ((t_up && ~t_up_prev) || (t_down && t_up_prev))...
                && ~firstAdjustment && nReversalsT >= nReversals(1)...
                && tStepModePos ~= length(stepModes) 
                switchTStepSize = true;
        end
        
        % If step sizes are locked to only be adjusted together, modify 
        % accordingly
        if pairStepSizes
            switchPStepSize = (switchPStepSize && switchTStepSize);
            switchTStepSize = switchPStepSize;
        end 
    end
    
    % ******** Process the keypresses *************
    % Match
    if match
        if ~silent
            Snd('Play',sin(0:5000));
        end
        % Calculate match. If this is an adjustment match, simply provide
        % the subject setting. For a forced-choice match, average the last
        % nReversals(2) light settings.
        if ~adjustment
            % Number of settings values to adjust
            [nSettings,~] = size(subjectSettings);
            nSettingsToAvg = min(nSettings,nReversals(2));
            
            % Average the selected number of values
            pMatch = mean(subjectSettings(end-(nSettingsToAvg-1):end,2));
            tMatch = mean(subjectSettings(end-(nSettingsToAvg-1):end,1));
            
            % Find the match positions
            pPos = zeros(1,nSettingsToAvg);
            tPos = zeros(1,nSettingsToAvg);
            for i = 0:(nSettingsToAvg-1)
                pPos(i+1) = find(p1Scales==subjectSettings(end-i,2));
                tPos(i+1) = find(testScales==subjectSettings(end-i,1));
            end
            primaryPos = mean(pPos);
            testPos = mean(tPos);
        else
            pMatch = p1Scales(primaryPos);
            tMatch = testScales(testPos);
        end
        fprintf('User found match at %g test, %g primary\n',...
            tMatch,pMatch);
        matches = [matches;[tMatch,pMatch]];
        matchPositions = [matchPositions;[testPos,primaryPos]];
        
        if plotResponses
            nAdjustments = length(subjectSettings(matchSettingInd:end,1));
            % Individual trajectory figure
            figure(nPlots);
            subplot(2,1,1)
            line5p = plot(nAdjustments+1,matches(end,2),'m* ',...
                nAdjustments+1,idealPRatio,'gs',...
                'MarkerSize',6,'LineWidth',1);
            if ~isempty(pStepSwitchInds)
                for i = 1:max(1,length(pStepSwitchInds))
                    line4p = line([pStepSwitchInds(i) pStepSwitchInds(i)],[0 1]);
                end
                legend([line5p(1:2)', line1p, line2p, line4p, line3p],...
                    'Subject Match','Nominal Match','Subject Settings',...
                    'Subject Responses','Step Size Changes','Subject Reversals',...
                    'Location','northeastoutside');
            else
                legend([line5p(1:2)', line1p, line2p, line3p],...
                    'Subject Match','Nominal Match','Subject Settings',...
                    'Subject Responses', 'Subject Reversals',...
                    'Location','northeastoutside');
            end 
                    
            
            subplot(2,1,2);
            line5t = plot(nAdjustments+1,matches(end,1),'m*',...
                nAdjustments+1,idealTestIntensity,'gs',...
                'MarkerSize',6,'LineWidth',1);
            if ~isempty(tStepSwitchInds)
                for i = 1:max(1,length(tStepSwitchInds))
                    line4t = line([tStepSwitchInds(i) tStepSwitchInds(i)],[0 1]);
                end
                legend([line5t(1:2)',line1t, line2t,line4t,line3t],...
                'Subject Match','Nominal Match','Subject Settings', ...
                'Subject Responses','Step Size Changes', 'Subject Reversals', ...
                'Location','northeastoutside');
            else
                legend([line5t(1:2)',line1t, line2t, line3t],...
                'Subject Match','Nominal Match','Subject Settings',...
                'Subject Responses', 'Subject Reversals',...
                'Location','northeastoutside');
            end
            
            % Joint trajectory figure
            figure(nPlots+1);
            hold on;
            p3 = plot(matches(end,2),matches(end,1),'r* ',...
                'MarkerSize',6,'LineWidth',1);
            p4 = plot(idealPRatio,idealTestIntensity,'gs',...
                'MarkerFaceColor','g');
            legend([p3 p4],'Subject Match','Nominal Match','Location','northeastoutside');
            
            % Prepare for next match
            if savePlots
                NicePlot.exportFigToPDF(fullfile(outputDir,...
                    [subjectID  '_' num2str(sessionNum) '_'...
                    num2str(size(matches,1)) '_plot']),plot1,300);
                NicePlot.exportFigToPDF(fullfile(outputDir,....
                    [subjectID '_' num2str(sessionNum)...
                    '_' num2str(size(matches,1)) '_jointPlot']),plot2,300);
            end
            matchSettingInd = length(subjectSettings(:,1))+1;
        end
        
        % Reset initial condition
        firstAdjustment = true;
    end
    
    % Quit
    if size(matches,1) == nObserverMatches    % Check if all matches have been made
        quit = true;
    end
    if quit
        if ~silent
            Snd('Play',sin(0:5000));
        end
        if plotResponses
            fprintf('User exited program\n');
        end
        save(fileLoc, 'matches','matchPositions',...
            'subjectSettings','p1','p2','test','calName',...
            'primarySpdsNominal','primarySpdsPredicted',...
            'testSpdsNominal','testSpdsPredicted',...
            'primaryStartStops','testStartStops',...
            'subjectID','allowedLambdaInds','allowedTIInds',...
            'sessionNum','annulusData','sInterval','lInterval',...
            'adjustmentLength','fieldSize','simNominalLights',...
            'nReversals','observer',...
            'p1Scale','p2Scale','testScale','lightFileName',...
            'adjustment','thresholdScaleFactor',...
            'nBelowThreshold','simObserver',...
            'p1Scales','testScales','observerParams',...
            'opponentParams','nObserverMatches','pIdealIndex',...
            'tIdealIndex','idealTestSpd','idealPrimarySpd',....
            'idealTestIntensity','idealPRatio','monochromatic','S');
        stillLooping = false;
        break;
    end
    
        % Switch p1 step size
    if switchPStepSize
        pStepModePos = pStepModePos+1;
        nReversalsP = 0;
        pStepSwitchInds = [pStepSwitchInds, displayLoopCounter];
        if pStepModePos > length(stepModes)
            pStepModePos = 1;
        end
        % Number of beeps indicates new step size position
        if ~silent
            for i = 1:pStepModePos
                Snd('Play',sin(0:5000));
            end
        end
        if plotResponses
            fprintf('User switched primary step size to %g\n',...
                (stepModes(pStepModePos)/(adjustmentLength-1)));
        end
    end
    
    % Switch test step size
    if switchTStepSize
        tStepModePos = tStepModePos+1;
        nReversalsT = 0;
        tStepSwitchInds = [tStepSwitchInds, displayLoopCounter];
        if tStepModePos > length(stepModes)
            tStepModePos = 1;
        end
        % Number of beeps indicates new step size position
        if ~silent
            for i = 1:tStepModePos
                Snd('Play',sin(0:5000));
            end
        end
        if plotResponses
            fprintf('User switched reference step size to %g\n',...
                (stepModes(tStepModePos)/(adjustmentLength-1)));
        end
    end
    
    % P1 up
    if p1_up
        primaryPos = primaryPos+stepModes(pStepModePos);
        if primaryPos > allowedLambdaInds(end)
            primaryPos = allowedLambdaInds(end);
            if ~silent
                Snd('Play',sin(0:5000));
            end
            if plotResponses
                fprintf('User reached upper primary limit\n');
            end
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                testScales(testPos),p1Scales(primaryPos));
        end
        % P1 down
    elseif p1_down
        primaryPos = primaryPos-stepModes(pStepModePos);
        if primaryPos < allowedLambdaInds(1)
            primaryPos = allowedLambdaInds(1);
            if ~silent
                Snd('Play',sin(0:5000));
            end
            if plotResponses
                fprintf('User reached lower primary limit\n');
            end
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                testScales(testPos),p1Scales(primaryPos));
        end
    end
    primaryRes = [primaryRes, double(p1_up)];
    
    % test up
    if t_up
        testPos = testPos+stepModes(tStepModePos);
        if testPos > allowedTIInds(end)
            testPos = allowedTIInds(end);
            if ~silent
                Snd('Play',sin(0:5000));
            end
            if plotResponses
                fprintf('User reached upper test limit\n');
            end
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                testScales(testPos),p1Scales(primaryPos));
        end
        % test down
    elseif t_down
        testPos = testPos-stepModes(tStepModePos);
        if testPos < allowedTIInds(1)
            testPos = allowedTIInds(1);
            if ~silent
                Snd('Play',sin(0:5000));
            end
            if plotResponses
                fprintf('User reached lower test limit\n');
            end
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                testScales(testPos),p1Scales(primaryPos));
        end
    end
    testRes = [testRes, double(t_up)];
    
    % Store data and set up for the next iteration
    p1_up_prev = p1_up;
    t_up_prev = t_up;
    subjectSettings = [subjectSettings;...
        [testScales(testPos),p1Scales(primaryPos)]];
    if (p1_up || p1_down || t_up || t_down)
        firstAdjustment = false;
    end
    
    % Edit plots
    if plotResponses
        figure(nPlots);
        subplot(2,1,1);
        line1p = plot(subjectSettings(matchSettingInd:end,2),'bo-','LineWidth',1);
        line2p = plot(primaryRes(matchSettingInd:end),...
            'r* ','MarkerSize',1);
        line3p = plot(pRevIndices(pRevIndices >= matchSettingInd)-matchSettingInd+1,...
            subjectSettings(pRevIndices(pRevIndices >= matchSettingInd),2),...
            'bo ', 'MarkerFaceColor','Blue');
        
        subplot(2,1,2);
        line1t = plot(subjectSettings(matchSettingInd:end,1),'bo-','LineWidth',1);
        line2t = plot(testRes(matchSettingInd:end),...
            'r* ','MarkerSize',1);
        line3t = plot(tRevIndices(tRevIndices >= matchSettingInd)-matchSettingInd+1,...
            subjectSettings(tRevIndices(tRevIndices >= matchSettingInd),1),...
            'bo ', 'MarkerFaceColor','Blue');
        
        figure(nPlots+1);
        plot(subjectSettings(matchSettingInd:end,2),...
            subjectSettings(matchSettingInd:end,1), 'bo-','LineWidth',1);
    end
end

%% Close up
if ~simObserver
    ol.setAll(false);
    if fieldSize > 2
        GLW_CloseAnnularStimulus();
    end
end

% Close extra plots
if plotResponses && simObserver
    close(figure(nPlots),figure(nPlots+1));
end
end