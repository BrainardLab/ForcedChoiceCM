function error = findMatchError(paramsVec,initialObs,testSpd,primarySpd)
% Computes the error associated with a given Rayleigh match 
%
% Syntax:
%   findMatchError(paramsVec,initialObs,testSpd,primarySpd)
%
% Description:
%    Given observer cone parameters, computes the "error" associated with a 
%    given pair of spectra which the observer identified as a Rayleigh 
%    match. Computes the observer's cone responses to both test and primary
%    lights, then converts this to opponent contrast. The error is
%    represented as the vector length of the luminance and RG opponent
%    contrast terms. 
%
%    This function is designed for use with fmincon or
%    similar parameter search functions. 
%
% Inputs:
%    paramsVec   -Vector of eight individual difference parameters. See 
%                 ObserverVecToParams for a full description  
%    initialObs  -Struct containing the initial settings for the observer.
%                 Not modified by the program, but some fields are used
%                 for reference. 
%    testSpd     -201x1 vector representation of the predicted spd for the
%                 chosen test light
%    primarySpd  -201x1 vector representation of the predicted spd for the
%                 chosen primary light
%
% Outputs:
%    error       -Vector length of the luminance and RG opponent contrast
%                 terms.
%
% Optional key-value pairs:
%    None 

% History:
%   06/12/20  dce       Wrote it.

% Generate a model observer with the given parameters 
params = [paramsVec initialObs.colorDiffParams.noiseSd]; % Append noise
observer = genRayleighObserver('age', initialObs.coneParams.ageYears,...
    'fieldSize', initialObs.coneParams.fieldSizeDegrees,'coneVec',params);

% Calculate cone responses for the given spectra 
test_LMS = observer.T_cones * testSpd; 
primary_LMS = observer.T_cones * primarySpd; 

% Calculate opponent contrast
opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
    test_LMS, primary_LMS);

% Take root mean square error of luminance and RG opponent components
error = norm(opponentContrast(1:2)); 
end 