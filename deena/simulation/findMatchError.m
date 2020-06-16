function error = findMatchError(paramsVec,initialObs,testSpds,primarySpds)
% Computes the error associated with a set of Rayleigh matches
%
% Syntax:
%   findMatchError(paramsVec,initialObs,testSpds,primarySpds)
%
% Description:
%    Given observer cone parameters, computes the "error" associated with
%    a set of pairs of spectra which the observer identified as Rayleigh
%    matches. For each pair, computes the observer's cone responses to both
%    test and primary lights, then converts this to opponent contrast. The
%    error for a given pair is represented as the vector length of the
%    luminance and RG opponent contrast terms. The overall error reported
%    is the root mean square of the individual error terms.
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
%    testSpd     -201 x n vector representation of the predicted spds for
%                 the chosen test light.
%    primarySpd  -201 x n vector representation of the predicted spds for
%                 the chosen primary light.
%
% Outputs:
%    error       -Root mean square of the vector lengths of luminance and
%                 RG opponent contrast terms.
%
% Optional key-value pairs:
%    None

% History:
%   06/12/20  dce       Wrote it.
%   06/15/20  dce       Modified to take in multiple spds

% Generate a model observer with the given parameters
params = [paramsVec 0]; % Append noise
observer = genRayleighObserver('age', initialObs.coneParams.ageYears,...
    'fieldSize', initialObs.coneParams.fieldSizeDegrees,'coneVec',params);

% How many match pairs do we have?
[~,nMatches] = size(testSpds);
pairError = zeros(1,nMatches);

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