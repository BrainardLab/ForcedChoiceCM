function [primaryCones, testCones] = OLGetConeEffects(fName, varargin)
% Calculate cone effects of subjects' OneLight Rayleigh matches.
%
% Syntax:
%   OLGetConeEffects(fName, varargin)
%
% Description
%    Takes in a file of user's Rayleigh matches on the OneLight (produced
%    by OLRayleighMatch or OLRayleighExp). For each match, calculates the
%    cone responses to the test and match lights. Can also run to plot test
%    match cone responses based on the radiometer playback of matches.
%
% Inputs:
%    fName         - character array of filename. Ends in .mat
%
% Outputs:
%    primaryCones  - 3xn matrix with cone responses for primary lights,
%                    where n is the number of matches. L cone responses are
%                    stored in the first row, M in the second, and S in the
%                    third.
%    testCones     - 3xn matrix with cone responses for test lights, where
%                    n is the number of matches. L cone responses are
%                    stored in the first row, M in the second, and S in the
%                    third.
%
% Optional key-value pairs:
%    'measured'  - logical indicating whether to calculate cone effects
%                  from radiometer measurements (true) or nominal spd data
%                  (false). Default is false.

% History
%    1/22/20   dce  - Modified program from OLTestConeEffects
%    4/5/20    dce  - Edited for style
% Example: OLGetConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/test/test_1.mat')

%% Parse input
close all;
p = inputParser;
p.addParameter('measured', false, @(x) (islogical(x)));
p.parse(varargin{:});
measured = p.Results.measured;

%% Load data file and radiometer file if appropriate
theData = load(fName);
if (isfield(theData,'matches') == 0 ||...
        isfield(theData,'matchPositions') == 0 ||...
        isfield(theData, 'cal') == 0 ||...
        isfield(theData, 'primarySpdsNominal') == 0 ||...
        isfield(theData,'primarySpdsPredicted') == 0 ||...
        isfield(theData, 'testSpdsNominal') == 0||...
        isfield(theData, 'testSpdsPredicted') == 0 ...
        || isfield(theData, 'foveal' == 0))
    error('Data file does not contain required variables');
end

if measured
    fprintf('\n******** Loading radiometer file ********\n');
    measFile = [fName(1:end - 4), '_meas', '.mat'];
    if ~exist(measFile, 'file')
        error('Radiometer measurements not available');
    end
    measData = load(measFile);
    if (isfield(measData,'measuredTestSpds') == 0 ||...
            isfield(measData,'measuredPrimarySpds') == 0)
        error('Radiometer file does not contain required variables');
    end
    fprintf('Radiometer measurements successfully loaded\n');
end

%% Initial Calculations
[nMatches, ~] = size(theData.matches); % Number of matches
wls = theData.cal.computed.pr650Wls;   % Wavelength range tested
inc = wls(2) - wls(1);                 % Increment of wavelengths

% Generate standard cone fundamentals for observer
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc,...
    'foveal', theData.foveal);

% Initialize arrays for storing data
primaryCones = zeros(3, nMatches);
testCones = zeros(3, nMatches);

%% Calculate effects of spectra on cones
% For each match, multiply the spd (predicted or measured) by the cone
% fundamentals
for i = 1:nMatches
    if ~measured
        testSpdPredicted = theData.testSpdsPredicted(:,theData.matchPositions(i,1));
        testCones(:,i) = T_cones * testSpdPredicted;
        
        primarySpdPredicted = theData.primarySpdsPredicted(:,theData.matchPositions(i,2));
        primaryCones(:,i) = T_cones * primarySpdPredicted;
    else
        testCones(:,i) = T_cones * measData.measuredTestSpds(:, i);
        primaryCones(:,i) = T_cones * measData.measuredPrimarySpds(:, i);
    end
end
end