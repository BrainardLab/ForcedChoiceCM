function comparisonLMS = OpponentContrastToLMS(colorDiffParams,referenceLMS,opponentContrast)
%
% Syntax:
%     comparisonLMS = OpponentContrastToLMS(colorDiffParams,referenceLMS,opponentContrast)
%
% Description:
%     Convert opponent contrast respresentation to LMS. This inverts
%     LMSToOpponentContrast.
%
% Inputs:
%     colorDiffParams        - Structure of parameters understood by LMSToOpponentContrast
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
% See also: LMSToOpponentContrast, ComputeMatchDiff
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
    opponentContrast = LMSToOpponentContrast(colorDiffParams,referenceLMS,comparisonLMS)
    checkLMS = OpponentContrastToLMS(colorDiffParams,referenceLMS,opponentContrast)
    if (max(abs(comparisonLMS - checkLMS) ./ comparisonLMS) > 1e-6)
        error('Routines do not self invert properly.');
    end
%}

switch (colorDiffParams.type)
    case 'opponentContrast'
        
        % Build up matrix
        M = diag([colorDiffParams.lumWeight colorDiffParams.rgWeight colorDiffParams.byWeight])*[ [colorDiffParams.LMRatio 1 0]/(colorDiffParams.LMRatio+1) ; ...
            [1 -1 0] ; ...
            [-0.5 -0.5 1] ];
        
        % Go to cone contrast
        coneContrast = inv(M)*opponentContrast;
        
        % And then LMS
        comparisonLMS = (coneContrast .* referenceLMS) + referenceLMS;
        
    otherwise
        error('Unknown type field passed in colorDiffParams structure');
end

