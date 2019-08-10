function [colorDiff,comparisonContrast] = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%
% Syntax:
%     [colorDiff,comparisonContrast] = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%
% Description:
%     Compute a single number color difference between two vectors
%     specified in cone coordinates. 
%
%     A vector length is taken in an opponent contrast representation,
%     which itself is computed by LMSToOpponentContrast.
%
% Inputs:
%     colorDiffParams        - Structure understood LMSToOpponentContrast
%     referenceLMS           - LMS coordinates of the reference with
%                              respect to which contrast is computed.
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
    colorDiff = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%}

% Get opponent representatoin
comparisonContrast = LMSToOpponentContrast(colorDiffParams,referenceLMS,comparisonLMS);

% Take appropriate weighted vector length
colorDiff = norm(comparisonContrast);
        
 
end

