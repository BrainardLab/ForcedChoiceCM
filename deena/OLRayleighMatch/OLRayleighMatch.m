function OLRayleighMatch(subjectID,sessionNum,varargin)
% Program to run a Rayleigh match experiment on the OneLight, or as a
% simulation
%
% Syntax:
%   OLRayleighMatch(subjectID,sessionNum)
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
%    judge the primary ratio and the test light intensity using the game
%    pad's directional pad. Judgements are made relative to the appearance
%    of the second light - right means the second light is comparatively
%    redder, left means the second light is comparatively greener, up means
%    the second light is comparatively brighter, and down means the second
%    light is comparatively dimmer. However, the program determines when to 
%    change the step size and when to declare a match. The forced-choice 
%    matches are an average of the last few settings, with the precise 
%    number specified by nReversals(2). 
%
%    In the "adjustment" procedure, the subject can use the top 'Y' button
%    to toggle primary step size, the left 'X' button to toggle test step
%    size, the right 'B' button to record a match, the 'Back' button to
%    show the ideal match, and the 'Start' button for switching light
%    order, in addition to the controls described above. Here, the
%    directional pad refers to desired adjustments, not appearance
%    judgements. Right makes the second light comparatively
%    redder, left makes the second light comparatively greener, up makes
%    the second light comparatively brighter, and down makes the second
%    light comparatively dimmer. In both procedures, the bottom 'A' button
%    can be used to force quit. The live experiment also includes an option
%    to run two interleaved Rayleigh matches at once - one with the primary
%    mixture displayed first, one with the reference light displayed first.
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
%    'testScale'       - Numerical scale factor for the test light, between
%                       0 and 1. Default is 0.1.
%    'isInterval'      - Interstimulus interval in s. Default is 0.1.
%    'itInterval'      - Intertrial interval in s. Default is 1.
%    'stimInterval'    - length of time that each stimulus is shown for, in
%                        s. Default is 0.5.
%    'fieldSize'      - Numeric field size, in degrees. Default is 2.
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
%    'nReversals'     - Number of reversals required before changing step 
%                       size, for simulated or forced choice observer. Entered
%                       as a 2-element vector - the first element is the
%                       number of reversals for intermediate step sizes,
%                       the second is the number needed for the smallest
%                       step size. Default is [1 4].
%    'adjustment'     - Logical indicating to use an adjustment method, not
%                       forced choice. For simulated observers, this means
%                       that matches are recorded when the cone excitation
%                       difference for the 2 spectra falls below a
%                       threshold. For live observers, this means that
%                       observers can freely change step sizes and make a
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
%     'coneNoise'            -Logical. If true, adds noise at the level of 
%                             the cones and not at the level of the opponent
%                             response. Default is false.
%    'nObserverMatches'     - Number of matches to simulate/run. Default is
%                             1.
%    'adjustmentLength'     - Number of possible steps available for
%                             adjusting the primary and test lights.
%                             Default is 201.
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
%    'testFirst'            -Logical. If true, displays the test light
%                            before the primary mixture, instead of the
%                            other way around. Default is false.
%    'pairStepSizes'        -Logical. If true, adjusts primary and test
%                            step sizes together instead of separately.
%                            Default is false.
%    'whiteScaleFactor'     -Scale factor for the neutral (white) light,
%                            between 0 and 1. Default is 0.001.
%    'interleaveStaircases' -Logical. If true, runs 2 interleaved
%                            staircases at once in a live experiment - one
%                            with primary first and one with reference 
%                            first. Note that this is separate from the 
%                            number of matches specified by
%                            nObserverMatches - that is the number of non-
%                            interleaved runs. Default is false. 
%    'outerFilename'        -Character vector name for an optional extra
%                            layer of hierarchy in the file saving tree
%                            (immediately before subject ID). If empty, the
%                            extra level of hierarchy is not included.
%                            Default is [].

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
%   02/05/21  dce       Fixed charcode system for using real OneLight
%   03/07/21  dce       Restructured to incorporate a forced-choice
%                       procedure and make judgements for both lights on
%                       each trial
%   03/08/21  dce       Added flicker back in while waiting for subject
%                       decisions
%   03/15/21  dce       Edited for style, added two new key-value pairs
%   03/25/21  dce       Changed keypress meanings, added white light
%   04/07/21  dce       Changed step size reduction rule
%   04/28/21  dce       Saved additional parameters
%   06/01/21  dce       Changed speech, saved match setting indices,
%                       added option to interleave staircases, edited for style
%   06/02/21  dce       Edited for style, got rid of monochromatic
%                       and nominal spd matching, fixed interleaving
%   06/04/21  dce       Step size adjustment, fixed plotting
%   06/04/21  dce       Added option to add noise to cone responses
%   06/21/21  dce       Fixed recording of reversal placement, changed
%                       match to be calculated based on last reversals 

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
p.addParameter('isInterval',0.1,@(x)(isnumeric(x)));
p.addParameter('itInterval',1,@(x)(isnumeric(x)));
p.addParameter('stimInterval',0.5,@(x)(isnumeric(x)));
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
p.addParameter('coneNoise',false,@(x)(islogical(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('lambdaRef',[],@(x)(isnumeric(x)));
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.addParameter('testFirst',false,@(x)(islogical(x)));
p.addParameter('pairStepSizes',false,@(x)(islogical(x)));
p.addParameter('whiteScaleFactor',0.001,@(x)(isnumeric(x)));
p.addParameter('outerFileName',[]);
p.addParameter('interleaveStaircases',false,@(x)(islogical(x)));
p.parse(varargin{:});

p1 = p.Results.p1;
p2 = p.Results.p2;
test = p.Results.test;
p1Scale = p.Results.p1Scale;
p2Scale = p.Results.p2Scale;
testScale = p.Results.testScale;
adjustmentLength = p.Results.adjustmentLength;
isInterval = p.Results.isInterval;
itInterval = p.Results.itInterval;
stimInterval = p.Results.stimInterval;
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
stimLimits = p.Results.stimLimits;
lambdaRef = p.Results.lambdaRef;
pairStepSizes = p.Results.pairStepSizes;
whiteScaleFactor = p.Results.whiteScaleFactor;
interleaveStaircases = p.Results.interleaveStaircases;
coneNoise = p.Results.coneNoise;

% Basic input error checking
if (length(nReversals)~=2)
    error('Reversal vector must be 2x1');
end
if (length(observerParams)~=8)
    error('Observer parameters must be entered as an 8-element vector');
end
if (length(opponentParams)~=4)
    error('Opponent contrast parameters must be entered as an 4-element vector');
end

%% Set up directory for saving results
% Create directory named subjectID for saving data, if it doesn't exist
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',p.Results.outerFileName,subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Create data file with name subject ID_SessionNumber. Throw error if a
% file already exists with that name.
fileLoc = fullfile(outputDir,[subjectID,'_',num2str(sessionNum),'.mat']);
if exist(fileLoc,'file')
    error('Specified output file %s already exists',[subjectID,'_',num2str(sessionNum),'.mat']);
end

%% Find light settings
% Use precomputed spectra, or compute if they do not exist
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
calName = lightSettings.cal.describe.calID;
primarySpds = lightSettings.primarySpdsPredicted;
testSpds = lightSettings.testSpdsPredicted;
primaryStartStops = lightSettings.primaryStartStops;
testStartStops = lightSettings.testStartStops;
p1Scales = lightSettings.p1Scales;
testScales = lightSettings.testScales;
S = lightSettings.cal.computed.pr650S;

whitePrimary = whiteScaleFactor * ones(size(lightSettings.cal.computed.pr650M,2),1);
[whiteStarts, whiteStops] = OLSettingsToStartsStops(lightSettings.cal,OLPrimaryToSettings(lightSettings.cal,whitePrimary));

% Find reference spd (for comparative contrasts), if desired
refSpd = [];
if ~isempty(lambdaRef)
    refSpd = primarySpds(:,p1Scales==lambdaRef);
    if isempty(refSpd)
        error('Provided reference lambda is not a selected primary mixture scalar');
    end
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

%% Initialize simulated observer, and find nominal match 
% If a cone vector was not specified, generates a standard observer
observer = genRayleighObserver('fieldSize',fieldSize,'age',age,...
    'coneVec',observerParams,'opponentParams',opponentParams,...
    'S',S);

% Find nominal match for observer. 
[~,~,tIdealIndex,pIdealIndex] =...
    searchPredictedRayleighMatch(testSpds,primarySpds,observer,'refSpd',refSpd);
idealTestIntensity = testScales(tIdealIndex);
idealPRatio = p1Scales(pIdealIndex);

%% Intialize OneLight and button box for live experiment 
if ~simObserver
    ol = OneLight();
    gamePad = GamePad();
end

%% Set up projector (if not making foveal or simulated matches)
annulusData = []; % Placeholder variable
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
% Possible step sizes
stepModes = [floor(adjustmentLength/5),floor(adjustmentLength/20),...
    floor(adjustmentLength/100),floor(adjustmentLength/200),...
    floor(adjustmentLength/400),floor(adjustmentLength/800),...
    floor(adjustmentLength/1600),floor(adjustmentLength/3200)];
stepModes = stepModes(stepModes ~= 0);

% Light display times 
waitTimes = [stimInterval isInterval stimInterval itInterval];

% Set up a data-storing structure. This holds data for a single Rayleigh
% match - for interleaved matching, use two copies of this struct.
matchData = struct();
matchData.ideal = false;              % Do not start with the ideal match 
matchData.displayLoopCounter = 0;     % Count number of iterations
matchData.p1_up_prev = true;          % Previous primary setting
matchData.t_up_prev = true;           % Previous ref setting
matchData.firstAdjustment = true;     % This is the first adjustment of the match
matchData.matches = [];               % Output array with subject matches
matchData.matchPositions = [];        % Positions of matches in the adjustment array
matchData.subjectSettings = [];       % Subject settings
matchData.primaryRes = [];            % Fill with 1's for every "up" trial, 0's for every "down" trial
matchData.testRes = [];               % Fill with 1's for every "up" trial, 0's for every "down" trial
matchData.pRevIndices = [];           % Indices where a reversal occurred
matchData.tRevIndices = [];           % Indices where a reversal occurred
matchData.pStepSwitchInds = [];       % Indices where primary step size was changed
matchData.tStepSwitchInds = [];       % Indices where ref step size was changed
matchData.matchSettingInds = [1];     % Indices where new match started

% Setup for if you use two interleaved staircases
if interleaveStaircases
    nTrialsPerLoop = 2;
    dataArr = {matchData matchData};
    dataArr{1}.nPlots = -3;
    dataArr{2}.nPlots = -1;
    staircaseTestFirst = [true false];   % Is the ref or the primary shown first?
    staircaseCompleted = [false false];  % Has the staircase been completed?
else
    nTrialsPerLoop = 1;
    dataArr = {matchData};
    dataArr{1}.nPlots = -1;
    staircaseTestFirst = p.Results.testFirst;
    staircaseCompleted = false;
end
% Track which staircase came first on a given trial. Enter 1 if the
% first staircase came first, 0 otherwise
superTrialOrderings = [];

%% Display loop
stillLooping = true;
while(stillLooping)
    % Choose how to order the interleaved trials in this supertrial
    if interleaveStaircases
        staircaseOrder = randperm(2);
    else 
        staircaseOrder = 1;
    end
    superTrialOrderings = [superTrialOrderings, staircaseOrder(1)];
    
    % Conduct a supertrial (one trial of each staircase)
    for tt = 1:nTrialsPerLoop
        % If all staircases are finished, the experiment is done. If
        % this staircase is completed, move to the next.
        if all(staircaseCompleted)
            stillLooping = false;
            break;
        elseif staircaseCompleted(staircaseOrder(tt))
            continue;
        else
            matchData = dataArr{staircaseOrder(tt)};
            testFirst = staircaseTestFirst(staircaseOrder(tt));
        end
        
        % Reset loop control parameters to initial states
        matchData.displayLoopCounter = matchData.displayLoopCounter + 1;
        match = false;
        quit = false;
        switchTStepSize = false;
        switchPStepSize = false;
        p1_up = false;
        p1_down = false;
        t_up = false;
        t_down = false;
        
        % Setup if it's a new match
        if matchData.firstAdjustment
            matchData.pStepModePos = 1;         % Start with the largest step sizes
            matchData.tStepModePos = 1;
            matchData.nReversalsP = 0;          % Initially, have made 0 reversals for each
            matchData.nReversalsT = 0;          % adjustment
            matchData.countBelowThreshold = 0;  % Count pairs below match threshold
            % Reset lights to a random position 
            matchData.primaryPos = randsample(allowedLambdaInds,1);
            matchData.testPos = randsample(allowedTIInds,1);
            matchData.subjectSettings = [matchData.subjectSettings;...
                [testScales(matchData.testPos),p1Scales(matchData.primaryPos)]];
            
            if plotResponses
                matchData.nPlots = matchData.nPlots+2*nTrialsPerLoop;
                % Plot with primary and test settings separately
                figure(matchData.nPlots);
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
                figure(matchData.nPlots+1);
                hold on;
                xlim([0 1]);
                ylim([0 1]);
                ylabel('Proportional Reference Intensity');
                xlabel('Primary Ratio (Proportion Red)');
                title2 = sprintf('Subject Settings Over Time, Staircase %g, Match %g',...
                    staircaseOrder(tt),ceil((matchData.nPlots+1-staircaseOrder(tt))/2));
                title(title2);
                plot(idealPRatio,idealTestIntensity,'gs',...
                    'MarkerFaceColor','g');
            end
        end
        % Define light indices
        if matchData.ideal
            pI = pIdealIndex;
            tI = tIdealIndex;
        else
            pI = matchData.primaryPos;
            tI = matchData.testPos;
        end
        
        if testFirst
            starts = {squeeze(testStartStops(tI,1,:))',whiteStarts,...
                squeeze(primaryStartStops(pI,1,:))',whiteStarts};
            stops = {squeeze(testStartStops(tI,2,:))',whiteStops, ...
                squeeze(primaryStartStops(pI,2,:))', whiteStops};
        else
            starts = {squeeze(primaryStartStops(pI,1,:))',whiteStarts,...
                squeeze(testStartStops(tI,1,:))',whiteStarts};
            stops = {squeeze(primaryStartStops(pI,2,:))',whiteStops,...
                squeeze(testStartStops(tI,2,:))',whiteStops};
        end
        
        % In a forced-choice live experiment, we show the specified lights then
        % pause for a decision period. In a simulated experiment, we prompt the
        % simulated observer for a decision.
        if (~simObserver)
            % For forced choice, we display the lights first before
            % prompting for the observer decision 
            if ~adjustment 
                for i = 1:4
                    nowTime = mglGetSecs;
                    lightInd = mod(i,4);
                    if lightInd == 0
                        lightInd = 4;
                    end
                    while(mglGetSecs < nowTime+waitTimes(lightInd)) 
                        ol.setMirrors(starts{lightInd},stops{lightInd});
                    end
                end
            end
            % Get the observer decision. If using the forced-choice method, we
            % get one R/G and one brightness decision. If using the adjustment
            % method, we wait for the first keypress.
            innerLoopCounter = 0;
            waitingForResponse = true;
            if ~adjustment         
                Speak('Redder or greener?'); % Prompt for first decision
            end
            while(waitingForResponse)
                %  Identify which light to display on this loop
                innerLoopCounter = innerLoopCounter+1;
                lightInd = mod(innerLoopCounter,4);
                if lightInd == 0
                    lightInd = 4;
                end
                nowTime = mglGetSecs;
                % Display light while waiting for observer decision
                while(mglGetSecs < nowTime+waitTimes(lightInd))
                    ol.setMirrors(starts{lightInd},stops{lightInd});
                    key = gamePad.getKeyEvent();
                    if (~isempty(key))
                        % Key interpretations - adjustment
                        if adjustment
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
                                    matchData.ideal = ~matchData.ideal;
                                    if ~silent
                                        Snd('Play',sin(0:5000));
                                    end
                                    if plotResponses
                                        if matchData.ideal
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
                                        fprintf('User switched order of primary and reference lights\n');
                                    end
                                case 'GP:North'
                                    if testFirst
                                        t_down = true;
                                    else
                                        t_up = true;
                                    end
                                case 'GP:South'
                                    if testFirst
                                        t_up = true;
                                    else
                                        t_down = true;
                                    end
                                case 'GP:East'
                                    if testFirst
                                        p1_up = true;
                                    else
                                        p1_down = true;
                                    end
                                case 'GP:West'
                                    if testFirst
                                        p1_down = true;
                                    else
                                        p1_up = true;
                                    end
                            end
                            waitingForResponse = false;
                            break;
                            
                        else % Key interpretations - FC red/green decision
                            switch(key.charCode)
                                case 'GP:East'
                                    if testFirst
                                        p1_down = true;
                                    else
                                        p1_up = true;
                                    end
                                    waitingForResponse = false;
                                    break;
                                case 'GP:West'
                                    if testFirst
                                        p1_up = true;
                                    else
                                        p1_down = true;
                                    end
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
            end
            if ~silent
                Snd('Play',sin(0:5000)/100);
            end
            
            % If doing an FC match, collect the light/dark decision
            if ~adjustment
                innerLoopCounter = 0;
                Speak('Brighter or dimmer?');
                waitingForResponse = true;
                while(waitingForResponse)
                    nowTime = mglGetSecs;
                    innerLoopCounter = innerLoopCounter+1;
                    lightInd = mod(innerLoopCounter,4);
                    if lightInd == 0
                        lightInd = 4;
                    end
                    while(mglGetSecs < nowTime+waitTimes(lightInd))  % Flicker while waiting for response
                        ol.setMirrors(starts{lightInd},stops{lightInd});
                        key = gamePad.getKeyEvent();
                        if (~isempty(key))
                            switch(key.charCode)
                                case 'GP:North'
                                    if testFirst
                                        t_up = true;
                                    else
                                        t_down = true;
                                    end
                                    waitingForResponse = false;
                                    break;
                                case 'GP:South'
                                    if testFirst
                                        t_down = true;
                                    else
                                        t_up = true;
                                    end
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
                if ~silent
                    Snd('Play',sin(0:5000)/100);
                end
            end
        else  % Prompt simulated observer for decisions
            [p1_up,t_up,isBelowThreshold] = ...
                observerRayleighDecision(observer,primarySpds(:,matchData.primaryPos),...
                testSpds(:,matchData.testPos),'thresholdScale',thresholdScaleFactor,...
                'noiseScale',noiseScaleFactor,'refSpd',refSpd,'coneNoise',...
                coneNoise);
            if ~p1_up
                p1_down = true;
            end
            if ~t_up
                t_down = true;
            end
            
            %  A couple methods specific for simulation with adjustment method.
            if adjustment && isBelowThreshold
                matchData.countBelowThreshold = matchData.countBelowThreshold+1;
            end
        end
        
        % Check for a primary reversal
        if ((p1_up && ~matchData.p1_up_prev) || (p1_down && matchData.p1_up_prev)) && ~matchData.firstAdjustment % Reversal
            matchData.nReversalsP = matchData.nReversalsP+1;  % Count the reversal
            matchData.pRevIndices = [matchData.pRevIndices,matchData.displayLoopCounter];
            if plotResponses
                fprintf('Primary reversal %g, step size %g, staircase %g\n',...
                    matchData.nReversalsP,(stepModes(matchData.pStepModePos)/...
                    (adjustmentLength-1)),staircaseOrder(tt));
            end
        end
        
        % Check for a ref reversal
        if ((t_up && ~matchData.t_up_prev) || (t_down && matchData.t_up_prev)) && ~matchData.firstAdjustment % Reversal
            matchData.nReversalsT = matchData.nReversalsT+1;  % Count the reversal
            matchData.tRevIndices = [matchData.tRevIndices,matchData.displayLoopCounter];
            if plotResponses
                fprintf('Test reversal %g, step size %g, staircase %g\n',...
                    matchData.nReversalsT,(stepModes(matchData.tStepModePos)/...
                    (adjustmentLength-1)),staircaseOrder(tt));
            end
        end
        
        % Process the observer responses, if not using the live adjustment
        % experiment
        if ~(adjustment && ~simObserver)
            % Check if we're at a match
            if (simObserver && adjustment && matchData.countBelowThreshold ==...
                    nBelowThreshold) || (~adjustment ...
                    && matchData.pStepModePos == length(stepModes)...
                    && matchData.tStepModePos == length(stepModes)...
                    && (matchData.nReversalsP >= nReversals(2))...
                    && (matchData.nReversalsT >= nReversals(2)))
                match = true;
            end
            % Check if primary step size needs to be adjusted
            if ((p1_up && ~matchData.p1_up_prev) || (p1_down && matchData.p1_up_prev))...
                    && ~matchData.firstAdjustment && matchData.nReversalsP >= nReversals(1)...
                    && matchData.pStepModePos ~= length(stepModes)
                switchPStepSize = true;
            end
            
            % Check if test step size needs to be adjusted
            if ((t_up && ~matchData.t_up_prev) || (t_down && matchData.t_up_prev))...
                    && ~matchData.firstAdjustment && matchData.nReversalsT >= nReversals(1)...
                    && matchData.tStepModePos ~= length(stepModes)
                switchTStepSize = true;
            end
            
            % If step sizes are locked to only be adjusted together, modify
            % accordingly. The current rule changes both step sizes if both
            % primary ratio and test intensity have undergone enough reversals
            % to reduce the step size, and one of the two is directed to adjust
            if pairStepSizes
                switchPStepSize = ((switchPStepSize || switchTStepSize)&&...
                    matchData.nReversalsP >= nReversals(1) && ...
                    matchData.nReversalsT >= nReversals(1));
                switchTStepSize = switchPStepSize;
            end
        end
        
        % ******** Process the keypresses *************
        % Switch p1 step size
        if switchPStepSize
            matchData.pStepModePos = matchData.pStepModePos+1;
            matchData.nReversalsP = 0;
            matchData.pStepSwitchInds = [matchData.pStepSwitchInds,matchData.displayLoopCounter];
            if matchData.pStepModePos > length(stepModes)
                matchData.pStepModePos = 1;
            end
            % Number of beeps indicates new step size position
            if ~silent && adjustment
                for i = 1:matchData.pStepModePos
                    Snd('Play',sin(0:5000));
                end
            end
            if plotResponses
                fprintf('User switched primary step size to %g for staircase %g \n',...
                    (stepModes(matchData.pStepModePos)/(adjustmentLength-1)), staircaseOrder(tt));
            end
        end
        
        % Switch test step size
        if switchTStepSize
            matchData.tStepModePos = matchData.tStepModePos+1;
            matchData.nReversalsT = 0;
            matchData.tStepSwitchInds = [matchData.tStepSwitchInds, matchData.displayLoopCounter];
            if matchData.tStepModePos > length(stepModes)
                matchData.tStepModePos = 1;
            end
            % Number of beeps indicates new step size position
            if ~silent && adjustment
                for i = 1:matchData.tStepModePos
                    Snd('Play',sin(0:5000));
                end
            end
            if plotResponses
                fprintf('User switched reference step size to %g for staircase %g\n',...
                    (stepModes(matchData.tStepModePos)/(adjustmentLength-1)),staircaseOrder(tt));
            end
        end
        
        % Announce step size change (works best for paired step sizes)
        if switchTStepSize
            announceStepSizeChange = true;
            for kk = 1:length(dataArr)
                if kk == staircaseOrder(tt)
                    continue;
                elseif dataArr{kk}.tStepModePos < matchData.tStepModePos
                    announceStepSizeChange = false;
                end
            end
            if announceStepSizeChange && ~silent
                stepString = sprintf('Step size %g of %g', matchData.tStepModePos, length(stepModes));
                Speak(stepString);
            end
        end
        
        % Adjust primary and reference (unless you're at the match point,
        % in which case further adjustment is irrelevant)
        matchData.testRes = [matchData.testRes, double(t_up)];
        matchData.primaryRes = [matchData.primaryRes, double(p1_up)];            
        if ~match
            % P1 up
            if p1_up
                matchData.primaryPos = matchData.primaryPos+stepModes(matchData.pStepModePos);
                if matchData.primaryPos > allowedLambdaInds(end)
                    matchData.primaryPos = allowedLambdaInds(end);
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    if plotResponses
                        fprintf('User reached upper primary limit for staicase %g\n',staircaseOrder(tt));
                    end
                end
                if plotResponses
                    fprintf('User pressed key for staircase %g. Test intensity = %g, red primary = %g \n',...
                        staircaseOrder(tt),testScales(matchData.testPos),p1Scales(matchData.primaryPos));
                end
                % P1 down
            elseif p1_down
                matchData.primaryPos = matchData.primaryPos-stepModes(matchData.pStepModePos);
                if matchData.primaryPos < allowedLambdaInds(1)
                    matchData.primaryPos = allowedLambdaInds(1);
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    if plotResponses
                        fprintf('User reached lower primary limit for staircase %g\n', staircaseOrder(tt));
                    end
                end
                if plotResponses
                    fprintf('User pressed key for staircase %g. Test intensity = %g, red primary = %g\n',...
                        staircaseOrder(tt),testScales(matchData.testPos),p1Scales(matchData.primaryPos));
                end
            end
            
            % test up
            if t_up
                matchData.testPos = matchData.testPos+stepModes(matchData.tStepModePos);
                if matchData.testPos > allowedTIInds(end)
                    matchData.testPos = allowedTIInds(end);
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    if plotResponses
                        fprintf('User reached upper test limit for staircase %g\n',staircaseOrder(tt));
                    end
                end
                if plotResponses
                    fprintf('User pressed key for staircase %g. Test intensity = %g, red primary = %g\n',...
                        staircaseOrder(tt),testScales(matchData.testPos),p1Scales(matchData.primaryPos));
                end
                % test down
            elseif t_down
                matchData.testPos = matchData.testPos-stepModes(matchData.tStepModePos);
                if matchData.testPos < allowedTIInds(1)
                    matchData.testPos = allowedTIInds(1);
                    if ~silent
                        Snd('Play',sin(0:5000));
                    end
                    if plotResponses
                        fprintf('User reached lower test limit for staircase %g\n',staircaseOrder(tt));
                    end
                end
                if plotResponses
                    fprintf('User pressed key for staircase %g. Test intensity = %g, red primary = %g\n',...
                        staircaseOrder(tt),testScales(matchData.testPos),p1Scales(matchData.primaryPos));
                end
            end
            
            % Store data and set up for the next iteration
            matchData.p1_up_prev = p1_up;
            matchData.t_up_prev = t_up;
            matchData.subjectSettings = [matchData.subjectSettings;...
                [testScales(matchData.testPos),p1Scales(matchData.primaryPos)]];
            if (p1_up || p1_down || t_up || t_down)
                matchData.firstAdjustment = false;
            end
        end
        
        % Edit plots
        if plotResponses
            figure(matchData.nPlots);
            subplot(2,1,1);
            line1p = plot(matchData.subjectSettings(matchData.matchSettingInds(end):end,2),...
                'bo-','LineWidth',1);
            line2p = plot(matchData.primaryRes(matchData.matchSettingInds(end):end),...
                'r* ','MarkerSize',1);
            line3p = plot(matchData.pRevIndices(matchData.pRevIndices >= matchData.matchSettingInds(end))...
                -matchData.matchSettingInds(end)+1,matchData.subjectSettings(...
                matchData.pRevIndices(matchData.pRevIndices >=...
                matchData.matchSettingInds(end)),2),'bo ', 'MarkerFaceColor','Blue');
            
            subplot(2,1,2);
            line1t = plot(matchData.subjectSettings(matchData.matchSettingInds(end):end,1),...
                'bo-','LineWidth',1);
            line2t = plot(matchData.testRes(matchData.matchSettingInds(end):end),...
                'r* ','MarkerSize',1);
            line3t = plot(matchData.tRevIndices(matchData.tRevIndices >= ...
                matchData.matchSettingInds(end))-matchData.matchSettingInds(end)+1,...
                matchData.subjectSettings(matchData.tRevIndices(...
                matchData.tRevIndices >= matchData.matchSettingInds(end)),1),...
                'bo ', 'MarkerFaceColor','Blue');
            
            figure(matchData.nPlots+1);
            plot(matchData.subjectSettings(matchData.matchSettingInds(end):end,2),...
                matchData.subjectSettings(matchData.matchSettingInds(end):end,1), 'bo-','LineWidth',1);
        end
        
        % Match
        if match
            if ~silent && adjustment
                Snd('Play',sin(0:5000));
            end
            % Calculate match. If this is an adjustment match, simply provide
            % the subject setting. For a forced-choice match, average the 
            % light settings from the last nReversals(2) reversals.
            if ~adjustment                
                % Average the selected number of reversals
                pMatch = mean(matchData.subjectSettings(...
                    matchData.pRevIndices(end-nReversals(2)+1:end),2));
                tMatch = mean(matchData.subjectSettings(...
                    matchData.tRevIndices(end-nReversals(2)+1:end),1));
                
                % Find the match positions
                pMatchPositions = zeros(1,nReversals(2));
                tMatchPositions = zeros(1,nReversals(2));
                for i = 1:nReversals(2)
                    pMatchPositions(i) = find(p1Scales==...
                        matchData.subjectSettings(matchData.pRevIndices(end-i+1),2));
                    tMatchPositions(i) = find(testScales==...
                        matchData.subjectSettings(matchData.tRevIndices(end-i+1),1));
                end
                pMatchPos = mean(pMatchPositions);
                tMatchPos = mean(tMatchPositions);
            else
                pMatch = p1Scales(matchData.primaryPos);
                tMatch = testScales(matchData.testPos);
                pMatchPos = matchData.primaryPos;
                tMatchPos = matchData.testPos;
            end
            fprintf('User found match at %g test, %g primary\n',...
                tMatch,pMatch);
            matchData.matches = [matchData.matches;[tMatch,pMatch]];
            matchData.matchPositions = [matchData.matchPositions;[tMatchPos,pMatchPos]];
            
            if plotResponses
                nAdjustments = length(matchData.subjectSettings(matchData.matchSettingInds(end):end,1));
                % Individual trajectory figure
                figure(matchData.nPlots);
                subplot(2,1,1)
                line5p = plot(nAdjustments+1,matchData.matches(end,2),'m* ',...
                    nAdjustments+1,idealPRatio,'gs',...
                    'MarkerSize',6,'LineWidth',1);
                trialPStepSwitchInds = matchData.pStepSwitchInds(matchData.pStepSwitchInds>=matchData.matchSettingInds(end));
                if ~isempty(trialPStepSwitchInds)
                    for i = 1:length(trialPStepSwitchInds)
                        line4p = line([trialPStepSwitchInds(i)-matchData.matchSettingInds(end)+1,...
                            trialPStepSwitchInds(i)-matchData.matchSettingInds(end)+1],[0 1]);
                    end
                    legend([line5p(1:2)', line1p, line2p, line3p, line4p],...
                        'Subject Match','Nominal Match','Subject Settings',...
                        'Subject Responses','Subject Reversals','Step Size Changes',...
                        'Location','northeastoutside');
                else
                    legend([line5p(1:2)', line1p, line2p, line3p],...
                        'Subject Match','Nominal Match','Subject Settings',...
                        'Subject Responses', 'Subject Reversals',...
                        'Location','northeastoutside');
                end
                
                subplot(2,1,2);
                line5t = plot(nAdjustments+1,matchData.matches(end,1),'m*',...
                    nAdjustments+1,idealTestIntensity,'gs',...
                    'MarkerSize',6,'LineWidth',1);
                trialTStepSwitchInds = matchData.tStepSwitchInds(matchData.tStepSwitchInds>=...
                    matchData.matchSettingInds(end));
                if ~isempty(trialTStepSwitchInds)
                    for i = 1:length(trialTStepSwitchInds)
                         line4t = line([trialTStepSwitchInds(i)-matchData.matchSettingInds(end)+1,...
                            trialTStepSwitchInds(i)-matchData.matchSettingInds(end)+1],[0 1]);
                    end
                    legend([line5t(1:2)',line1t, line2t,line3t,line4t],...
                        'Subject Match','Nominal Match','Subject Settings', ...
                        'Subject Responses','Subject Reversals','Step Size Changes', ...
                        'Location','northeastoutside');
                else
                    legend([line5t(1:2)',line1t, line2t, line3t],...
                        'Subject Match','Nominal Match','Subject Settings',...
                        'Subject Responses', 'Subject Reversals',...
                        'Location','northeastoutside');
                end
                
                % Joint trajectory figure
                figure(matchData.nPlots+1);
                hold on;
                p3 = plot(matchData.matches(end,2),matchData.matches(end,1),'r* ',...
                    'MarkerSize',6,'LineWidth',1);
                p4 = plot(idealPRatio,idealTestIntensity,'gs',...
                    'MarkerFaceColor','g');
                legend([p3 p4],'Subject Match','Nominal Match','Location','northeastoutside');
                
                % Prepare for next match
                if savePlots
                    plot1 = figure(matchData.nPlots);
                    plot2 = figure(matchData.nPlots+1);
                    NicePlot.exportFigToPDF(fullfile(outputDir,...
                        [subjectID  '_' num2str(sessionNum) '_'...
                        num2str(size(matchData.matches,1)) '_' ...
                        num2str(staircaseOrder(tt)) '_plot']),plot1,300);
                    NicePlot.exportFigToPDF(fullfile(outputDir,....
                        [subjectID '_' num2str(sessionNum)...
                        '_' num2str(size(matchData.matches,1)) '_'...
                        num2str(staircaseOrder(tt)) '_jointPlot']),plot2,300);
                end
                matchData.matchSettingInds(end+1) = length(matchData.subjectSettings(:,1))+1;
            end
            
            % Reset initial condition
            matchData.firstAdjustment = true;
            if size(matchData.matches,1) == nObserverMatches    % Check if all matches have been made
                quit = true;
            end
        end
        dataArr{staircaseOrder(tt)} = matchData; % Update data before ending iteration

        % Quit
        if quit
            if plotResponses
                fprintf('User completed Rayleigh match %g\n',staircaseOrder(tt));
            end
            staircaseCompleted(staircaseOrder(tt)) = true;
            save(fileLoc,'dataArr','age','p1','p2','test','calName',...
                'primarySpds','testSpds','primaryStartStops','testStartStops',...
                'subjectID','allowedLambdaInds','allowedTIInds',...
                'sessionNum','annulusData','isInterval','itInterval','stimInterval',...
                'adjustmentLength','fieldSize','nReversals','observer',...
                'whiteScaleFactor','p1Scale','p2Scale','testScale',...
                'lightFileName','adjustment','thresholdScaleFactor',...
                'nBelowThreshold','simObserver','p1Scales','testScales',...
                'observerParams','opponentParams','nObserverMatches',...
                'pIdealIndex','tIdealIndex','idealTestIntensity','idealPRatio','S',...
                'pairStepSizes','whiteStarts','whiteStops','resetAnnulus',...
                'silent','staircaseTestFirst','lambdaRef','stimLimits','noiseScaleFactor',...
                'interleaveStaircases','superTrialOrderings','coneNoise');
        end
    end
end

%% Close up (live experiment)
if ~simObserver
    ol.setAll(false);
    if fieldSize > 2
        GLW_CloseAnnularStimulus();
    end
end
fprintf('User exited program\n');
if ~silent
    Snd('Play',sin(0:5000));
end
end