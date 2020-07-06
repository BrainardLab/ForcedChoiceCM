function error = findMatchError(paramsVec,initialObs,testSpds,primarySpds,varargin)
% Computes the error associated with a set of Rayleigh matches
%
% Syntax:
%   findMatchError(paramsVec,initialObs,testSpds,primarySpds)
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

% History:
%   06/12/20  dce       Wrote it.
%   06/15/20  dce       Modified to take in multiple spds

% Parse input 
p = inputParser;
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Generate a model observer with the given parameters
params = [paramsVec 0];          % Append noise
observer = genRayleighObserver('age', initialObs.coneParams.ageYears,...
    'fieldSize', initialObs.coneParams.fieldSizeDegrees,...
    'coneVec',params,'S',p.Results.S);

[spdLength,nMatches] = size(testSpds);   
% Throw error if matrix sizes do not match 
if length(SToWls(p.Results.S)) ~= spdLength
    error('Observer wavelength sampling and spd length do not match'); 
end 
pairError = zeros(1,nMatches);   % Array for storing error of each pair 

% Calculate opponent contrast for each match 
for i = 1:nMatches
    % Calculate cone responses for the given spectra
    test_LMS = observer.T_cones * testSpds(:,i);
    primary_LMS = observer.T_cones * primarySpds(:,i);
    
    % Calculate opponent contrast
    opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
        test_LMS, primary_LMS);
    
    % Find vector length of luminance and RG opponent components
    pairError(i) = norm(opponentContrast(1:2));
end

% Report the root mean squared error 
error = sqrt(mean(pairError.^2)); 
end