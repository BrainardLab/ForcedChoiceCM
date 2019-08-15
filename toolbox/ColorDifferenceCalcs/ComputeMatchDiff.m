function [colorDiff,comparisonContrast] = ComputeMatchDiff(colorDiffParams,M,adaptationLMS,referenceLMS,comparisonLMS)
%
% Syntax:
%     [colorDiff,comparisonContrast] = ComputeMatchDiff(colorDiffParams,M,adaptationLMS,referenceLMS,comparisonLMS)
%
% Description:
%     Compute a single number color difference between two vectors
%     specified in cone coordinates. 
%
%     A vector length of the difference between reference and comparison
%     is taken in an opponent contrast representation. The opponent
%     representation is defined by LMSToOpponentContrast.
%
% Inputs:
%     colorDiffParams        - Structure understood LMSToOpponentContrast
%     M                      - Matrix that goes from LMS contrast to opponent contrast
%     adpatationLMS          - LMS coordinates of the adapting feield
%     referenceLMS           - LMS coordinates of the reference.
%     comparisonLMS          - LMS coordinates of the stimulus whose
%                              contrast is computed.
%
%
% Outputs:
%   colorDiff                - Single number color difference measure.
%   comparisonContrast       - Contrast representation of the comparison
%
% Optional key/value pairs:
%   None:
%
% See also: LMSToOpponentContrast, OpponentContrastToLMS
%

% History:
%   08/09/19  dhb   Wrote it.

% Examples:
%{
    colorDiffParams.type = 'opponentContrast';
    colorDiffParams.LMRatio = 2;
    colorDiffParams.lumWeight = 4;
    colorDiffParams.rgWeight = 2;
    colorDiffParams.byWeight = 0.5;
    referenceLMS = [1 1 1]';
    comparisonLMS = [2 0.5 1.5]';
    M = GetOpponentContrastMatrix(colorDiffParams);
    colorDiff = ComputeMatchDiff(colorDiffParams,M,referenceLMS,referenceLMS,comparisonLMS)
%}

% Get opponent representations
referenceContrast = LMSToOpponentContrast(M,adaptationLMS,referenceLMS);
comparisonContrast = LMSToOpponentContrast(M,adaptationLMS,comparisonLMS);

% Take appropriate weighted vector length
colorDiff = norm(comparisonContrast-referenceContrast);
         
end

