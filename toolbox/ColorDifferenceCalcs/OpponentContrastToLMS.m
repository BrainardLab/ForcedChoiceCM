function comparisonLMS = OpponentContrastToLMS(M,referenceLMS,opponentContrast)
%
% Syntax:
%     comparisonLMS = OpponentContrastToLMS(M,referenceLMS,opponentContrast)
%
% Description:
%     Convert opponent contrast respresentation to LMS. This inverts
%     LMSToOpponentContrast.
%
% Inputs:
%     M                      - Matrix that goes from LMS contrast to
%                              opponent contrast.  See
%                              GetOpponentContrastMatrix.
%     referenceLMS           - LMS coordinates of the reference with
%                              respect to which contrast is computed.
%     opponentContrast       - Contrast representation to be converted.
%
% Outputs:
%     comparisonLMS          - Returned LMS representation.
%
% Optional key/value pairs:
%   None:
%
% See also: GetOpponentContrastMatrix, LMSToOpponentContrast, ComputeMatchDiff
%

% History:
%   08/09/19  dhb   Wrote it.

% Examples:
%{
    colorDiffParams.type = 'opponentContrast';
    colorDiffParams.LMRatio = 2;
    colorDiffParams.lumWeight = 1;
    colorDiffParams.rgWeight = 3;
    colorDiffParams.byWeight = 1.5;
    referenceLMS = [1 1 1]';
    comparisonLMS = [2 0.5 1.5]';
    M = GetOpponentContrastMatrix(colorDiffParams);
    opponentContrast = LMSToOpponentContrast(M,referenceLMS,comparisonLMS)
    checkLMS = OpponentContrastToLMS(M,referenceLMS,opponentContrast)
    if (max(abs(comparisonLMS - checkLMS) ./ comparisonLMS) > 1e-6)
        error('Routines do not self invert properly.');
    end
%}
        
% Go to cone contrast
coneContrast = inv(M)*opponentContrast;

% And then LMS
comparisonLMS = (coneContrast .* referenceLMS) + referenceLMS;

