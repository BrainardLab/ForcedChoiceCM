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
%                    Default is true.
%    'measLastOnly'- Logical. If true, measures only the last setting, not
%                    the last few reversals. Default is false.

% History
%    dce    xx/xx/19  - Wrote it
%    dce    3/29/20   - Edited for style, added ongoing saving
%    dce    2/13/21   - Restructured to measure matches from multiple files
%                       at once.
%    dce    2/25/21   - Added measurement of dark spd.
%    dce    6/01/21   - Restructure to explicitly measure the number of
%                       reversals needed to record a match
%    dce    6/02/21   - Edited to reflect changes in structure of
%                       OLRayleighMatch data files
%    dce    6/24/21   - Edited to calculate matches based on last
%                       reversals, not last settings
%    dce    6/28/21   - Added option to measure only the last setting, not
%                       the last few settings. This should be used in cases 
%                       of adjustment matching, but not forced choice.
%    dce    7/22/21   - Added saving of unaveraged measured spectra

%% Parse input
p = inputParser;
p.addParameter('measWhite', true, @(x) (islogical(x)));
p.addParameter('measLastOnly', false, @(x) (islogical(x)));
p.parse(varargin{:});

%% Set up output file for saving results
% Check a few options for what the output directory could be-start with a
% more specific subfolder, if not look for a more general subfolder
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjID,[subjID '_' num2str(sessionNum)]);
if ~exist(outputDir,'dir')
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'matchFiles',subjID);
    if ~exist(outputDir,'dir')       
        error('Match file not found for specified subject');
    end
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
    % Load data file
    theData = load(matchFiles{kk});   
    
    % Check if we are interleaving matches. If so, change array size
    % (We assume that either all data files contain interleaved matches, or
    % none do)
    if kk==1 && length(theData.dataArr)==2
        measuredPrimarySpdsUnavg = cell(2,length(matchFiles));
        measuredRefSpdsUnavg = cell(2,length(matchFiles));
    else
        measuredPrimarySpdsUnavg = cell(1,length(matchFiles));
        measuredRefSpdsUnavg = cell(1,length(matchFiles));
    end
    
    % Display and measure white light, if desired
    if p.Results.measWhite
        fprintf('Measuring white light...\n');
        ol.setMirrors(squeeze(theData.whiteStarts)', squeeze(theData.whiteStops)');
        pause(0.1);
        measuredWhite = [measuredWhite;spectroRadiometerOBJ.measure];
        fprintf('Successfully measured white light. Begin measuring matches...\n');
        save(outputFile,'measuredRefSpds', 'measuredPrimarySpds', 'measuredWhite');
    else
        measuredWhite = [measuredWhite, 0]; % Dummy variable for saving file
    end
    
    % Loop through matches within the file; display and measure each one
    for tt = 1:length(theData.dataArr) % Number of interleaved staircases
        [nMatches, ~] = size(theData.dataArr{tt}.matchPositions);
        measuredPrimarySpdsUnavg{tt,kk} = cell(1,nMatches);
        measuredRefSpdsUnavg{tt,kk} = cell(1,nMatches);
        for i = 1:nMatches % Number of matches within each staircase
            if p.Results.measLastOnly
                if i==nMatches
                    matchSettings = theData.dataArr{tt}.subjectSettings(end,:);
                else 
                    matchInd = theData.dataArr{tt}.matchSettingInds(i+1)-1;
                    matchSettings = theData.dataArr{tt}.subjectSettings(matchInd,:);
                end 
                nSettingsToMeasure = 1;
                pValsToMeasure = matchSettings(2);
                refValsToMeasure = matchSettings(1);
            else
                % Find settings and reversal indices we are using
                if i==nMatches  % Last match in the file
                    matchInds = theData.dataArr{tt}.matchSettingInds(i):...
                        size(theData.dataArr{tt}.subjectSettings,1);
                    pRevSettings = theData.dataArr{tt}.pRevIndices(...
                        theData.dataArr{tt}.pRevIndices>=theData.dataArr{tt}.matchSettingInds(i));
                    refRevSettings = theData.dataArr{tt}.tRevIndices(...
                        theData.dataArr{tt}.tRevIndices>=theData.dataArr{tt}.matchSettingInds(i));
                else
                    matchInds = theData.dataArr{tt}.matchSettingInds(i):...
                        theData.dataArr{tt}.matchSettingInds(i+1)-1;
                    pRevSettings = theData.dataArr{tt}.pRevIndices(...
                        (theData.dataArr{tt}.pRevIndices>=theData.dataArr{tt}.matchSettingInds(i))...
                        & (theData.dataArr{tt}.pRevIndices<theData.dataArr{tt}.matchSettingInds(i+1)));
                    refRevSettings = theData.dataArr{tt}.tRevIndices(...
                        (theData.dataArr{tt}.tRevIndices>=theData.dataArr{tt}.matchSettingInds(i))...
                        & (theData.dataArr{tt}.tRevIndices<theData.dataArr{tt}.matchSettingInds(i+1)));
                end
                matchSettings = theData.dataArr{tt}.subjectSettings(matchInds,:);
                
                % Find which indices we are measuring. By default, the
                % last nReversals(2) reversals are measured and averaged. If
                % fewer reversals are available, the last nReversals(2)
                % settings are measured. If fewer settings are available, then
                % all available settings are measured and averaged.
                if length(pRevSettings) < theData.nReversals(2) || length(refRevSettings) < theData.nReversals(2)
                    nSettingsToMeasure = min(theData.nReversals(2),size(matchSettings,1));
                    pValsToMeasure = matchSettings(end-nSettingsToMeasure+1:end,2);
                    refValsToMeasure = matchSettings(end-nSettingsToMeasure+1:end,1);
                else
                    nSettingsToMeasure = theData.nReversals(2);
                    pValsToMeasure = matchSettings(pRevSettings(end-nSettingsToMeasure+1:end),2);
                    refValsToMeasure = matchSettings(refRevSettings(end-nSettingsToMeasure+1:end),1);
                end
            end 
            
            primaryMeas = [];
            refMeas = [];
            for j = 1:nSettingsToMeasure
                trialPInd = find(theData.p1Scales==pValsToMeasure(j));
                trialRefInd = find(theData.testScales==refValsToMeasure(j));
                
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
            
            measuredPrimarySpdsUnavg{tt,kk}{i} = primaryMeas;
            measuredRefSpdsUnavg{tt,kk}{i} = refMeas;
            measuredPrimarySpds = [measuredPrimarySpds; mean(primaryMeas,1)];
            measuredRefSpds = [measuredRefSpds; mean(refMeas,1)];
            save(outputFile, 'measuredRefSpds', 'measuredPrimarySpds', 'measuredWhite',...
                'measuredPrimarySpdsUnavg','measuredRefSpdsUnavg');
        end
    end    
    fprintf('File %g of %g complete\n',kk,length(matchFiles));
end

% Measure dark spd
ol.setAll(false);
pause(0.1);
measuredDarkSpd = spectroRadiometerOBJ.measure;
save(outputFile, 'measuredRefSpds', 'measuredPrimarySpds', 'measuredWhite',...
    'measuredDarkSpd','measuredPrimarySpdsUnavg','measuredRefSpdsUnavg');

%% Close devices
spectroRadiometerOBJ.shutDown;
ol.setAll(false);
end