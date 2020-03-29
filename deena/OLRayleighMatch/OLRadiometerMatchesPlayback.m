function OLRadiometerMatchesPlayback(matchFile, varargin)
% Play back OneLight Rayleigh matches and measure with radiometer
%
% Syntax:
%   OLRayleighMatch(matchFile)
%
% Description
%    Takes in a file of user's Rayleigh matches on the OneLight. Plays
%    these back and records using the radiometer.
%
% Inputs:
%    matchfile  - character array of filename. Ends in .mat
%
% Outputs:
%    saves a file named 'matchfile_meas.mat' in the same directory as
%    matchfile.
%
% Optional key-value pairs:
%    'measWhite'   - Logical indicating whether to measure white light.
%                    Default is false.
%    'measIdeal'   - Logical indicating whether to measure the ideal match.
%                    Default is false.

% Example: OLRadiometerMatchesPlayback('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatches/test/test_1.mat')
% History
%    dce    xx/xx/19  - Wrote it
%    dce    3/29/20   - Edited for style, added ongoing saving 

%% Parse input
p = inputParser;
p.addParameter('measWhite', false, @(x) (islogical(x)));
p.addParameter('measIdeal', false, @(x) (islogical(x)));
p.parse(varargin{:});

%% Load file and check it contains required variables
theData = load(matchFile);
if(isfield(theData,'matchPositions') == 0 || ...
        isfield(theData,'matches') == 0 || ...
        isfield(theData,'primarySpdsPredicted') == 0 || ...
        isfield(theData,'primaryStartStops') == 0 ||...
        isfield(theData,'testStartStops') == 0 || ...
        isfield(theData,'foveal') == 0)
    error('Data file does not contain required variables');
elseif p.Results.measWhite &&...
        (isfield(theData,'whiteStarts') == 0 || isfield(theData, 'whiteStops' == 0))
    error('Data file does not contain variables for white light');
end

%% Set up radiometer
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

%% Initial setup 
ol = OneLight();
fName = [matchFile(1:end - 4), '_meas', '.mat']; % Filename to save 

[numMatches, ~] = size(theData.matches);
[spdLength, ~] = size(theData.primarySpdsPredicted);
measuredPrimarySpds = zeros(spdLength, numMatches);
measuredTestSpds = zeros(spdLength, numMatches);

%% Display and measure white light
if p.Results.measWhite
    fprintf('Measuring white light...\n');
    ol.setMirrors(squeeze(theData.whiteStarts)', squeeze(theData.whiteStops)');
    pause(0.1);
    measuredWhite = spectroRadiometerOBJ.measure;
    fprintf('Successfully measured white light. Begin measuring matches...\n');
    save(fName, 'measuredWhite');
else
    measuredWhite = 0; % Dummy variable for saving file
end

%% Measure ideal spectrum
if p.Results.measIdeal
    % Ideal indices--need to check on OneLight
    if theData.foveal
        pIdealIndex = 198;
        tIdealIndex = 12;
    else
        pIdealIndex = 199;
        tIdealIndex = 10;
    end
    fprintf('Measuring ideal light...\n');
    ol.setMirrors(squeeze(theData.primaryStartStops(pIdealIndex,1,:))',...
        squeeze(theData.primaryStartStops(pIdealIndex,2,:))');
    pause(0.1);
    measuredIdealPrimary = spectroRadiometerOBJ.measure;
    save(fName, 'measuredWhite', 'measuredIdealPrimary');
        
    ol.setMirrors(squeeze(theData.testStartStops(tIdealIndex,1,:))',...
        squeeze(theData.testStartStops(tIdealIndex,2,:))');
    measuredIdealTest = spectroRadiometerOBJ.measure;
    save(fName, 'measuredWhite', 'measuredIdealPrimary', 'measuredIdealTest');
else 
    % Set variables to zero so program will save 
    measuredIdealPrimary = 0; 
    measuredIdealTest = 0; 
end

%% Loop through matches; display and measure each one
for i = 1:numMatches
    % Display primary on OL, measure with radiometer
    ol.setMirrors(squeeze(theData.primaryStartStops(theData.matchPositions(i,2),1,:))',...
        squeeze(theData.primaryStartStops(theData.matchPositions(i,2),2,:))');
    pause(0.1);
    primaryMeas = spectroRadiometerOBJ.measure;
    measuredPrimarySpds(:,i) = primaryMeas;
    save(fName, 'measuredTestSpds', 'measuredPrimarySpds', 'measuredWhite',...
    'measuredIdealPrimary', 'measuredIdealTest');
    
    % Display test on OL, measure with radiometer
    ol.setMirrors(squeeze(theData.testStartStops(theData.matchPositions(i,1),1,:))',...
        squeeze(theData.testStartStops(theData.matchPositions(i,1),2,:))');
    pause(0.1);
    testMeas = spectroRadiometerOBJ.measure;
    measuredTestSpds(:,i) = testMeas;
    save(fName, 'measuredTestSpds', 'measuredPrimarySpds', 'measuredWhite',...
    'measuredIdealPrimary', 'measuredIdealTest');
    
    fprintf('Match %g complete\n', i);
end

% Close devices
spectroRadiometerOBJ.shutDown;
ol.setAll(false);
end