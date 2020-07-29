function err = opponentParamsErr(opponentParams,refLMS,shiftedLMS,noiseSD)
% Computes an error metric for the difference between opponent contrast and
% noise standard deviation.
%
% Syntax:
%   opponentParamsErr(opponentParams,refLMS,shiftedLMS,noiseSD)
%
% Description:
%    Takes in cone responses to a reference light (refLMS), as well as cone
%    responses to a set of lights shifted around the reference
%    (shiftedLMS). Using the provided opponent contrast weighting
%    metrics, computes the opponent contrast of the reference light and
%    each shifted light, then compares the opponent contrast norm to the
%    standard deviation. The returned error value is the rmse of the
%    difference between the opponent contrast norm and the noise standard 
%    deviation for each shifted light.
%
%    This function is designed to be used with parameter search routines
%    such as fmincon. It can be used to find the parameters that bring the 
%    opponent contrast scaling in line with CIELAB coordinates, if the 
%    shiftedLMS values are chosen to all have the same CIELAB distance from 
%    the reference.
%
% Inputs:
%    opponentParams -3x1 vector of opponent contrast parameters.
%                    (1) is the luminance weight, (2) is the rg weight, and
%                    (3) is the by weight.
%    refLMS         -3x1 vector of reference cone responses. 
%    shiftedLMS     -3xn matrix of cone responses for a series of test
%                    lights varied around the reference light.
%    noiseSD        -Measure of noise standard deviation (often set to 1)
%
% Outputs:
%    err           -Root mean square error of the difference between the
%                   noise standard deviation and the norm of the opponent 
%                   contrast terms.
%
% Optional key-value pairs:
%    None          

% History:
%   07/29/20  dce       Wrote it.

% Generate color difference params struct 
colorDiffParams = struct();
colorDiffParams.type = 'opponentContrast';
colorDiffParams.LMRatio = 2;
colorDiffParams.lumWeight = opponentParams(1);
colorDiffParams.rgWeight = opponentParams(2);
colorDiffParams.byWeight = opponentParams(3);
colorDiffParams.M = GetOpponentContrastMatrix(colorDiffParams);

[~,nPoints] = size(shiftedLMS); % Number of shifted points provided
indErr = zeros(nPoints,1);         % Vector for storing error information 
for i = 1:nPoints
    % Convert LMS values to opponent contrast 
    opponentContrast = LMSToOpponentContrast(colorDiffParams,refLMS,shiftedLMS);
    indErr(i) = norm(opponentContrast)-noiseSD; 
end
% Find rms error of the original error values 
err = sqrt(mean(indErr.^2));