% t_FitAsanoToRayleigh
%
% Description:
%    Demonstration/test tutorial to check/develop fitting methods to
%    Rayleigh match data.
%

% History:
%    07/29/21  dhb  Started.

%% Clear
clear; close all;

%% Load and homogenize data
whichData = 'simulatedNoise';
theDir = fileparts(mfilename('fullpath'));
switch (whichData)
    case 'simulatedNoise'
        dataFilename = 'SimulatedWithNoise_spds.mat';
        theData = load(fullfile(theDir,'sampleData',dataFilename));
        S = [380 2 201]; wls = SToWls(S);
        fieldSize = theData.fieldSize;
        age = theData.age;
        [nWls,nReps,nRefs] = size(theData.primarySpdsByWl);
        refSpds = reshape(theData.refSpdsByWl,S(3),nReps*nRefs);
        primarySpds = reshape(theData.primarySpdsByWl,S(3),nReps*nRefs);
    case 'MELC_0004'
        dataFilename = 'MELC_0004_spds.mat';
        theData = load(fullfile(theDir,'sampleData',dataFilename));
        S = [380 2 201]; wls = SToWls(S);
        fieldSize = 10;
        age = 32;
        [nWls,nReps,nRefs] = size(theData.measPrimarySpdsByWl);
        refSpds = reshape(theData.measRefSpdsByWl,S(3),nReps*nRefs);
        primarySpds = reshape(theData.measPrimarySpdsByWl,S(3),nReps*nRefs);
    otherwise
        error('What is this data file of which you speak?')
end

%% Call Deena fit routine
% findObserverParameters(testSpds,primarySpds,varargin)
params = findObserverParameters(refSpds,primarySpds,'S0',true,'dmac0',true,'sdLambdaMax',4, ...
    'age',age,'fieldSize',fieldSize, ...
    'minimizeConeErr',false,'matchErrorFun',@findMatchErrorNew);

params


function error = findMatchErrorNew(coneParamsVec,initialObs,testSpds,...
    primarySpds,varargin)
% Computes the error associated with a set of Rayleigh matches
%
% Syntax:
%   findMatchErrorNew(paramsVec,initialObs,testSpds,primarySpds)
%
% Description:
%    Given observer cone parameters, computes the error associated with
%    a set of pairs of spectra which the observer identified as Rayleigh
%    matches. For each pair, computes the observer's cone responses to both
%    test and primary lights, then converts this to opponent contrast. The
%    error for a given pair is represented as the vector length of the
%    luminance and RG opponent contrast terms. The overall error is the 
%    root mean square of the individual error terms.
%
%    This function is designed for use with fmincon or similar parameter
%    search functions.
%
% Inputs:
%    paramsVec   -Vector of eight individual difference parameters. See
%                 ObserverVecToParams for a full description
%    initialObs  -Struct containing the initial settings for the observer.
%                 Not modified by the program, but some fields are used
%                 for reference.
%    testSpd     -Vector representation of the predicted spds for
%                 the chosen test light.
%    primarySpd  -Vector representation of the predicted spds for
%                 the chosen primary light.
%
% Outputs:
%    error       -Root mean square of the vector lengths of luminance and
%                 RG opponent contrast terms.
%
% Optional key-value pairs:
%    S           -Wavelength sampling for cone calculations, in the 
%                 form [start increment numTerms]. Default is [380 2 201] 
%                 (OneLight convention)  
%    errScalar   -Integer for scaling the match error, to improve search.
%                 Default is 100.
%    findConeErr -Logical. If true, calculates cone exictation error
%                 instead of opponent contrast difference. Default 
%                 is false.

% History:
%   06/12/20  dce       Wrote it.
%   06/15/20  dce       Modified to take in multiple spds
%   05/09/21  dce       Added option to find error based on cone excitation
%                       difference instead of opponent contrast.

% Parse input 
p = inputParser;
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.addParameter('findConeErr',false,@(x)(islogical(x)));
p.parse(varargin{:});

[spdLength,nMatches] = size(testSpds);   
% Throw error if matrix sizes do not match 
if length(SToWls(p.Results.S)) ~= spdLength
    error('Observer wavelength sampling and spd length do not match'); 
end 
pairError = zeros(1,nMatches);   % Array for storing error of each pair

% Find opponent parameters
CDParams = initialObs.colorDiffParams;
opponentParams = [CDParams.lumWeight CDParams.rgWeight CDParams.byWeight...
    CDParams.noiseSd];

% Generate an observer
observer = genRayleighObserver('age', initialObs.coneParams.ageYears,...
    'fieldSize', initialObs.coneParams.fieldSizeDegrees,...
    'coneVec',coneParamsVec,'opponentParams',opponentParams,'S',p.Results.S); 

% Calculate error metric for each match 
for i = 1:nMatches
    % Calculate cone responses for the given spectra
    test_LMS = observer.T_cones * testSpds(:,i);
    primary_LMS = observer.T_cones * primarySpds(:,i);
    
    if p.Results.findConeErr
        % Calculate cone contrast difference
        error = 0.5*(test_LMS - primary_LMS) ./ (test_LMS + primary_LMS);
    else
        % Calculate opponent contrast
        opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
            test_LMS, primary_LMS);
        error = opponentContrast;
    end

    % Find vector length of error (excluding S components)
    pairError(i) = norm(error(1:2));
end

% Report the root mean square error (scaled up to improve searching)
error = sqrt(mean(pairError.^2))*p.Results.errScalar; 
end