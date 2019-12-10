% Program to play back OneLight Rayleigh matches and measure with
% radiometer
function OLRadiometerMatchesPlayback(matchFile)
% Syntax:
%   OLRayleighMatch(matchFile)
%
% Description
%    Takes in a file of user's Rayleigh matches on the OneLight. Plays
%    these back and records using the radiometer 
%
% Inputs:
%    matchfile  - character array of filename. Ends in .mat
%
% Outputs:
%    saves a file named 'matchfile_meas.mat' in the same directory as
%    matchfile. 
%
% Optional key-value pairs:
%    none 

% Example: OLRadiometerMatchesPlayback('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatches/test/test_1.mat')

%% Load file and check it contains required variables
theData = load(matchFile);
if (isfield(theData,'p1') == 0 || isfield(theData,'p2') == 0  ||...
        isfield(theData,'test') == 0 || isfield(theData,'matches') == 0 ...
        || isfield(theData,'matchPositions') == 0 || isfield(theData, 'cal') == 0 ||...
        isfield(theData, 'primarySpdsNominal' == 0) || isfield(theData,...
        'primarySpdsPredicted' == 0) || isfield(theData,...
        'testSpdsNominal' == 0)|| isfield(theData, 'testSpdsPredicted' == 0) ...
        || isfield(theData,'primaryStartStops') == 0 || isfield(theData,'testStartStops') == 0)
    error('Data file does not contain required variables');
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

%% Set up initial files
ol = OneLight();
[numMatches, ~] = size(theData.matches); 
[spdLength, ~] = size(theData.primarySpdsPredicted); 
measuredPrimarySpds = zeros(spdLength, numMatches); 
measuredTestSpds = zeros(spdLength, numMatches); 

%% Loop through matches 
 
for i = 1:numMatches
    % Display primary on OL, measure with radiometer
    ol.setMirrors(squeeze(theData.primaryStartStops(theData.matchPositions(i,2),1,:))',...
            squeeze(theData.primaryStartStops(theData.matchPositions(i,2),2,:))');
    pause(0.1);
    primaryMeas = spectroRadiometerOBJ.measure;
    measuredPrimarySpds(:,i) = primaryMeas;  
    
    % Display test on OL, measure with radiometer
    ol.setMirrors(squeeze(theData.testStartStops(theData.matchPositions(i,1),1,:))',...
            squeeze(theData.testStartStops(theData.matchPositions(i,1),2,:))');
    pause(0.1);
    testMeas = spectroRadiometerOBJ.measure;
    measuredTestSpds(:,i) = testMeas; 
end

% Save data
fName = [matchFile(1:end - 4), '_meas', '.mat']; 
save(fName, 'measuredTestSpds', 'measuredPrimarySpds');
spectroRadiometerOBJ.shutDown;
disp('Successfully played back matches'); 
end 