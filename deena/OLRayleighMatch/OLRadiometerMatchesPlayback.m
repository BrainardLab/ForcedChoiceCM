function OLRadiometerMatchesPlayback(subjID,sessionNum,matchFiles,varargin)
% Play back OneLight Rayleigh matches and measure with radiometer
%
% Syntax:
%   OLRadiometerMatchesPlayback(subjID,sessionNum,matchFiles)
%
% Description
%    Takes in a cell array of files storing data on a user's Rayleigh 
%    matches on the OneLight. Plays these back and records using the 
%    radiometer. Saves results as a single file.
%
% Inputs:
%    subjID      - subject ID, entered as a character vector
%    sessionNum  - Integer session number
%    matchfiles  - cell array where each cell is a filename ending in .mat
%
% Outputs:
%    Saves a file named 'subjID_meas.mat' in a subject-specific directory.
%
% Optional key-value pairs:
%    'measWhite'   - Logical indicating whether to measure white light.
%                    Default is false.

% History
%    dce    xx/xx/19  - Wrote it
%    dce    3/29/20   - Edited for style, added ongoing saving
%    dce    2/13/21   - Restructured to measure matches from multiple files
%                       at once.
%    dce    2/25/21   - Added measurement of dark spd. 
%    dce    6/01/21   - Restructure to explicitly measure the number of
%                       reversals needed to record a match

%% Parse input
p = inputParser;
p.addParameter('measWhite', false, @(x) (islogical(x)));
p.parse(varargin{:});

%% Set up output file for saving results
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjID,[subjID '_' num2str(sessionNum)]);
if ~exist(outputDir,'dir')
    error('Match data not found for specified subject');
end
outputFile = fullfile(outputDir,[subjID '_' num2str(sessionNum) '_meas.mat']);

% Set up arrays for storing radiometer data 
measuredPrimarySpds = [];
measuredRefSpds = [];
measuredWhite = [];

%% Set up OneLight and radiometer
ol = OneLight();

spectroRadiometerOBJ  = PR670dev(...
    'verbosity',        1, ...       % 1 -> minimum verbosity
    'devicePortString', [] ...       % empty -> automatic port detection)
    );

% Set options Options available for PR670:
spectroRadiometerOBJ.setOptions(...
    'verbosity',        1, ...
    'syncMode',         'OFF', ...      % choose from 'OFF', 'AUTO', [20 400];
    'cyclesToAverage',  1, ...          % choose any integer in range [1 99]
    'sensitivityMode',  'STANDARD', ... % choose between 'STANDARD' and 'EXTENDED'.  'STANDARD': (exposure range: 6 - 6,000 msec, 'EXTENDED': exposure range: 6 - 30,000 msec
    'exposureTime',     'ADAPTIVE', ... % choose between 'ADAPTIVE' (for adaptive exposure), or a value in the range [6 6000] for 'STANDARD' sensitivity mode, or a value in the range [6 30000] for the 'EXTENDED' sensitivity mode
    'apertureSize',     '1 DEG' ...   % choose between '1 DEG', '1/2 DEG', '1/4 DEG', '1/8 DEG'
    );

%% Loop through match files. For each one, display and measure lights 
for kk = 1:length(matchFiles)
    
    % Load data file and check if it contains required fields
    theData = load(matchFiles{kk});
    if(isfield(theData,'matchPositions') == 0 || ...
            isfield(theData,'matches') == 0 || ...
            isfield(theData,'nReversals') == 0 || ...
            isfield(theData,'subjectSettings') == 0 || ...
            isfield(theData,'matchSettingInds') == 0 || ...
            isfield(theData,'primarySpdsPredicted') == 0 || ...
            isfield(theData,'primaryStartStops') == 0 ||...
            isfield(theData,'testSpdsPredicted') == 0 || ...
            isfield(theData,'testStartStops') == 0)
        error('Data file does not contain required variables');
    elseif p.Results.measWhite &&...
            (isfield(theData,'whiteStarts') == 0 || isfield(theData, 'whiteStops' == 0))
        error('Data file does not contain variables for white light');
    end
    
    % Display and measure white light, if desired 
    if p.Results.measWhite
        fprintf('Measuring white light...\n');
        ol.setMirrors(squeeze(theData.whiteStarts)', squeeze(theData.whiteStops)');
        pause(0.1);
        measuredWhite = [measuredWhite;spectroRadiometerOBJ.measure];
        fprintf('Successfully measured white light. Begin measuring matches...\n');
        save(outputFile,'measuredTestSpds', 'measuredPrimarySpds', 'measuredWhite');
    else
        measuredWhite = [measuredWhite, 0]; % Dummy variable for saving file
    end
    
    % Loop through matches within the file; display and measure each one
    [nMatches, ~] = size(theData.matchPositions);
    for i = 1:nMatches
        % Find settings indices we are using
        if i==nMatches
            matchInds = theData.matchSettingInds(i):size(theData.subjectSettings,1);
        else 
             matchInds = theData.matchSettingInds(i):theData.matchSettingInds(i+1)-1;
        end 
        matchSettings = theData.subjectSettings(matchInds);
        
        % Number of settings to average for a given match 
        nSettingsToAvg = min(size(matchSettings,1),theData.nReversals);
        
        primaryMeas = [];
        refMeas = [];
        for j = 0:(nSettingsToAvg-1)
            trialPInd = find(p1Scales==matchSettings(end-j,2));
            trialRefInd = find(testScales==matchSettings(end-j,1));
            
            % Display primary on OL, measure with radiometer
            ol.setMirrors(squeeze(theData.primaryStartStops(trialPInd,1,:))',...
                squeeze(theData.primaryStartStops(trialPInd,2,:))');
            pause(0.1);
            primaryMeas = [primaryMeas; spectroRadiometerOBJ.measure];
            
            % Display reference on OL, measure with radiometer
            ol.setMirrors(squeeze(theData.testStartStops(trialRefInd,1,:))',...
                squeeze(theData.testStartStops(trialRefInd,2,:))');
            pause(0.1);
            refMeas = [refMeas; spectroRadiometerOBJ.measure];
        end
        
        measuredPrimarySpds = [measuredPrimarySpds; mean(primaryMeas,1)];
        measuredRefSpds = [measuredRefSpds,mean(refMeas,1)];
        save(outputFile, 'measuredRefSpds', 'measuredPrimarySpds', 'measuredWhite');
    end
    fprintf('File %g of %g complete\n',kk,length(matchFiles));
end

% Measure dark spd
ol.setAll(false);
pause(0.1);
measuredDarkSpd = spectroRadiometerOBJ.measure;
save(outputFile, 'measuredRefSpds', 'measuredPrimarySpds', 'measuredWhite',...
    'measuredDarkSpd');

%% Close devices
spectroRadiometerOBJ.shutDown;
ol.setAll(false);
end