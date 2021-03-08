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
%    best match for this observer. When this option is used, the experimenter
%    can choose between two decision rules - a "forced choice" procedure that
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
%                       is 600.
%    'p1Scale'        - Numerical scale factor for the first primary light,
%                       between 0 and 1. Default is 1.
%    'p2Scale'        - Numerical scale factor for the second primary
%                       light, between 0 and 1. Default is 0.02.
%    'testScale'      - Numerical scale factor for the test light, between
%                       0 and 1. Default is 0.1.
%    'sInterval'      - length of time in s that the short stimulus is
%                       displayed for. Default is 0.5.
%    'lInterval'      - length of time in s that the long stimulus is
%                       displayed for. Default is 1.
%    'foveal'         - logical indicating whether we are making foveal
%                       matches, in which case the annulus is not turned
%                       on. Default is true.
%    'age'            - Integer age for observer. Default is 32.
%    'resetAnnulus'   - logical indicating to run a script that lets the
%                       experimenter reset the annulus. Default is false.
%    'plotResponses'  - Logical indicating to make plots of the
%                       experimental timecourse and print output. Default
%                       is true.
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
%                             Default is 3201.
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
p.addParameter('adjustmentLength',3201,@(x)(isnumeric(x)));
p.addParameter('sInterval',0.5,@(x)(isnumeric(x)));
p.addParameter('lInterval',1,@(x)(isnumeric(x)));
p.addParameter('foveal',true,@(x)(islogical(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('resetAnnulus',false,@(x)(islogical(x)));
p.addParameter('plotResponses',true,@(x) (islogical(x)));
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
foveal = p.Results.foveal;
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
simNominalLights = p.Results.simNominalLights;
monochromatic = p.Results.monochromatic;
monochromaticS = p.Results.monochromaticS;
stimLimits = p.Results.stimLimits;
lambdaRef = p.Results.lambdaRef;

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

if foveal
    fieldSize = 2;
else
    fieldSize = 10;
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
    if isempty(allowedTIInds) || length(allowedTIInds) <3
        [~,minInd] = min(abs(testScales - stimLimits(3)));
        [~,maxInd] = min(abs(p1Scales - stimLimits(4)));
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
% Compute the nominal match
idealForStandardObs = false;
if idealForStandardObs
    [idealTestSpd,idealPrimarySpd,tIdealIndex,pIdealIndex] =...
        searchPredictedRayleighmatch(testSpds,primarySpds,stdObs,'refSpd',refSpd);
else
    [idealTestSpd,idealPrimarySpd,tIdealIndex,pIdealIndex] =...
        searchPredictedRayleighMatch(testSpds,primarySpds,observer,'refSpd',refSpd);
end
idealTestIntensity = testScales(tIdealIndex);
idealPRatio = p1Scales(pIdealIndex);

%% Intialize OneLight and button box/keypresses
if ~simObserver
    ol = OneLight();
    gamePad = GamePad();
end

%% Set up projector (if not making foveal or simulated matches)
annulusData = []; % Empty placeholder variable
% Check if an annulus file exists. If it does, can choose whether to use
% the existing file or reset the annulus.
if ~simObserver && ~foveal
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

% Initial display settings
primaryPos = allowedLambdaInds(1); % Start with first allowed position in primary array
if monochromatic && primaryPos==1
    primaryPos = 2;
end
testPos = allowedTIInds(1);        % Start with first allowed position in test array
if monochromatic && testPos==1
    testPos = 2;
end
p1_up_prev = true;    % Indicate that the program directed p1 and t to
t_up_prev = true;     % increase on prior trials.

nReversalsP = 0;      % Initially, have made 0 reversals for each
nReversalsT = 0;      % adjustment

pStepModePos = 1;            % Start with the largest step sizes
tStepModePos = 1;

countBelowThreshold = 0;    % Count pairs below match threshold
rev = false;                % Start with forward, not reverse order
ideal = false;              % Do not start with the ideal match
stillLooping = true;        % Start looping
firstAdjustment = true;     % This is the first adjustment of the match

% Data-storing arrays
matches = [];               % Output array with subject matches
matchPositions = [];        % Positions of matches in the adjustment array
subjectSettings = [testScales(testPos) p1Scales(primaryPos)];

% Optional plotting setup
if plotResponses
    nPlots = 1;           % Counter for current figure index
    matchSettingInd = 1;  % SubjectSettings position for current match start
end

%% Display loop
% Loop through primary and test light until the user presses a key or the
% simulated observer decides to simulate a keypress
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
    
    % Define light indices
    if ideal
        pI = pIdealIndex;
        tI = tIdealIndex;
    else
        pI = primaryPos;
        tI = testPos;
    end
    
    % In a forced-choice live experiment, we show the specified lights then
    % pause for a decision period. In a simulated experiment, we prompt the
    % simulated observer for a decision.
    if ~simObserver
        % Display the two lights, then go to a dark screen
        if rev
            ol.setMirrors(squeeze(testStartStops(tI,1,:))',...
                squeeze(testStartStops(tI,2,:))');
            pause(sInterval);
            ol.setMirrors(squeeze(primaryStartStops(pI,1,:))',...
                squeeze(primaryStartStops(pI,2,:))');
        else
            ol.setMirrors(squeeze(primaryStartStops(pI,1,:))',...
                squeeze(primaryStartStops(pI,2,:))');
            pause(sInterval);
            ol.setMirrors(squeeze(testStartStops(tI,1,:))',...
                squeeze(testStartStops(tI,2,:))');
        end
        pause(sInterval);
        ol.setAll(false);
        
        % Get the observer decision. If using the forced-choice method, we
        % get one R/G and one brightness decision. If using the adjustment
        % method, we wait for the first keypress.
        waitingForResponse = true;
        if adjustment
            while(waitingForResponse)
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
                        case 'GP:Back'  % Change to another key
                            rev = ~rev;
                            if ~silent
                                % One beep = primary first, two beeps = test first
                                Snd('Play',sin(0:5000)/ 20);
                                if rev
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
                end
            end
            pause(lInterval);
            
        else   % Forced choice case
            % Prompt for R/G decision
            Speak('Redder or Greener?');
            while(waitingForResponse)
                key = gamePad.getKeyEvent();
                if (~isempty(key))
                    switch(key.charCode)
                        case 'GP:East'
                            p1_up = true;
                            waitingForResponse = false;
                        case 'GP:West'
                            p1_down = true;
                            waitingForResponse = false;
                        otherwise
                            Speak('Invalid key');
                    end
                end
            end
            pause(lInterval);
            
            % Prompt for ref decision
            Speak('Brighter or Darker?');
            waitingForResponse = true;
            while(waitingForResponse)
                key = gamePad.getKeyEvent();
                if (~isempty(key))
                    switch(key.charCode)
                        case 'GP:North'
                            t_up = true;
                            waitingForResponse = false;
                        case 'GP:South'
                            t_down = true;
                            waitingForResponse = false;
                        otherwise
                            Speak('Invalid key');
                    end
                end
            end
            pause(lInterval);
        end
    else  % Prompt simulated observer for decisions
        [pJudge,tJudge,isBelowThreshold] = ...
            observerRayleighDecision(observer,primarySpds(:,primaryPos),...
            testSpds(:,testPos),'thresholdScale',thresholdScaleFactor,...
            'noiseScale',noiseScaleFactor,'refSpd',refSpd);
        if pJudge
            p1_up = true;
        else
            p1_down = true;
        end
        if tJudge
            t_up = true;
        else 
            t_down = true;
        end
        
        %  A couple simulation - specific methods
        % If we are making threshold matches, record whether this
        % combination of lights is below the threshold
        if adjustment && isBelowThreshold
            countBelowThreshold = countBelowThreshold+1;
            if plotResponses
                fprintf('Below threshold: %g\n',countBelowThreshold);
            end
        end
        if adjustment && (nReversalsP >=...
                10*nReversals(2)) &&(nReversalsT >= 10*nReversals(2))
            % Quit if we're stuck in an infinite loop
            quit = true;
            fprintf('Terminated due to infinite loop\n');
        end
    end
    
    % Make new plots if it's the start of a match
    if plotResponses && firstAdjustment
        % Plot with primary and test settings separately
        figure(nPlots);
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
        figure(nPlots+1);
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
            p1_up_prev = logical(round(rand()));
            t_up_prev = logical(round(rand()));
            nReversalsP = 0;
            nReversalsT = 0;
            countBelowThreshold = 0;
        end
        
        % Check if primary step size needs to be adjusted
        if ((p1_up && ~p1_up_prev) || (p1_down && p1_up_prev)) && ~firstAdjustment % Reversal 
            if nReversalsP == nReversals(1) && pStepModePos ~= stepModes(end) % Lower step size
                switchPStepSize = true;
            else % Count the reversal
                nReversalsP = nReversalsP+1;
                if plotResponses
                    fprintf('Primary reversal %g, step size %g\n',...
                        nReversalsP,(stepModes(pStepModePos)/...
                        (adjustmentLength-1)));
                end
            end
        end
        p1_up_prev = p1_up;
        
        % Check if reference step size needs to be adjusted
        if ((t_up && ~t_up_prev) || (t_down && t_up_prev)) && ~firstAdjustment% Reversal 
            if nReversalsT == nReversals(1) && tStepModePos ~= stepModes(end) % Lower step size
                switchTStepSize = true;
            else % Count the reversal
                nReversalsT = nReversalsT+1;
                if plotResponses
                    fprintf('Reference reversal %g, step size %g\n',...
                        nReversalsT,(stepModes(tStepModePos)/...
                        (adjustmentLength-1)));
                end
            end
        end
        t_up_prev = t_up;
    end 
          
    % ******** Process the keypresses *************
    % Match
    if match
        if ~silent
            Snd('Play',sin(0:5000));
        end
        % Calculate match. If this is an adjustment match, simply provide 
        % the subject setting. For a forced-choice match, average the last 
        % four light settings.
        if ~adjustment
            % Number of settings values to adjust
            [nSettings,~] = size(subjectSettings);
            nSettingsToAvg = min(nSettings,4);
            
            % Average the selected number of values 
            pMatch = mean(subjectSettings(end-(nSettingsToAvg-1):end,2));
            tMatch = mean(subjectSettings(end-(nSettingsToAvg-1):end,1));
            
            % Find the match positions
            pPos = zeros(1,nSettingsToAvg);
            tPos = zeros(1,nSettingsToAvg);
            for i = 0:(nSettingsToAvg-1)
                pPos(i) = find(p1Scales==subjectSettings(end-i,2));
                tPos(i) = find(testScales==subjectSettings(end-i,1));
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
            p1 = plot(nAdjustments+1,matches(end,2),'r* ',...
                nAdjustments+1,idealPRatio,'gs',...
                'MarkerSize',6,'LineWidth',1);
            legend(p1,'Subject Match','Nominal Match');
            subplot(2,1,2);
            p2 = plot(nAdjustments+1,matches(end,1),'r*',...
                nAdjustments+1,idealTestIntensity,'gs',...
                'MarkerSize',6,'LineWidth',1);
            legend(p2,'Subject Match','Nominal Match');
            
            % Joint trajectory figure
            figure(nPlots+1);
            hold on;
            p3 = plot(matches(end,2),matches(end,1),'r* ',...
                'MarkerSize',6,'LineWidth',1);
            p4 = plot(idealPRatio,idealTestIntensity,'gs',...
                'MarkerFaceColor','g');
            legend([p3 p4],'Subject Match','Nominal Match');
            
            % Prepare for next match
            nPlots = nPlots+2;
            matchSettingInd = length(subjectSettings(:,1))+1;
        end
        
        % Reset initial condition, and reset lights randomly
        pStepModePos = 1;
        tStepModePos = 1;
        firstAdjustment = true;
        primaryPos = randsample(allowedLambdaInds,1);
        testPos = randsample(allowedTIInds,1);
        if monochromatic 
            if primaryPos==1
                primaryPos = 2;
            end
            if testPos==1
                testPos = 2;
            end
        end
        subjectSettings = [subjectSettings;...
            [testScales(testPos),p1Scales(primaryPos)]];
    end
    
    % Quit
    [row,~] = size(matches);
    if row == nObserverMatches    % Check if all matches have been made
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
            'adjustmentLength','foveal','simNominalLights',...
            'nReversals','observer',...
            'p1Scale','p2Scale','testScale','lightFileName',...
            'adjustment','thresholdScaleFactor',...
            'nBelowThreshold','simObserver',...
            'p1Scales','testScales','observerParams',...
            'opponentParams','nObserverMatches','pIdealIndex',...
            'tIdealIndex','idealForStandardObs','idealTestSpd',...
            'idealPrimarySpd','idealTestIntensity','idealPRatio',...
            'monochromatic','S');
        stillLooping = false;
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
        if ~silent
            Snd('Play',sin(0:5000)/100);
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g \n',...
                testScales(testPos),p1Scales(primaryPos));
        end
        subjectSettings = [subjectSettings;...
            [testScales(testPos),p1Scales(primaryPos)]];
        firstAdjustment = false;
    end
    
    % P1 down
    if p1_down
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
        if ~silent
            Snd('Play',sin(0:5000)/100);
        end
        subjectSettings = [subjectSettings;...
            [testScales(testPos),p1Scales(primaryPos)]];
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                testScales(testPos),p1Scales(primaryPos));
        end
        firstAdjustment = false;
    end
    
    % Ref up
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
        if ~silent
            Snd('Play',sin(0:5000)/100);
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                testScales(testPos),p1Scales(primaryPos));
        end
        subjectSettings = [subjectSettings;...
            [testScales(testPos),p1Scales(primaryPos)]];
        firstAdjustment = false;
    end
    
    % Ref down
    if t_down
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
        if ~silent
            Snd('Play',sin(0:5000)/100);
        end
        if plotResponses
            fprintf('User pressed key. Test intensity = %g, red primary = %g\n',...
                testScales(testPos),p1Scales(primaryPos));
        end
        subjectSettings = [subjectSettings;...
            [testScales(testPos),p1Scales(primaryPos)]];
        firstAdjustment = false;
    end
    
    % Switch p1 step size
    if switchPStepSize
        pStepModePos = pStepModePos+1;
        nReversalsP = 0;
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
    
    % Switch ref step size
    if switchTStepSize
        tStepModePos = tStepModePos+1;
        nReversalsT = 0;
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
    
    % Edit plots if an adjustment was made
    if plotResponses && (p1_up || p1_down || t_up || t_down)
        figure(nPlots);
        subplot(2,1,1);
        plot(subjectSettings(matchSettingInd:end,2),'bo-','LineWidth',1);
        subplot(2,1,2);
        plot(subjectSettings(matchSettingInd:end,1),'bo-','LineWidth',1);
        figure(nPlots+1);
        plot(subjectSettings(matchSettingInd:end,2),...
            subjectSettings(matchSettingInd:end,1), 'bo-','LineWidth',1);
    end 
end 
    
%% Close up
if ~simObserver
    ol.setAll(false);
    if ~foveal
        GLW_CloseAnnularStimulus();
    end
end
% Close extra plots
if plotResponses && simObserver
    close(figure(nPlots),figure(nPlots+1));
end
end