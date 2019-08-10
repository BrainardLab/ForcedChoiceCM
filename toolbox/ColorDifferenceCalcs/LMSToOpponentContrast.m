function opponentContrast = LMSToOpponentContrast(colorDiffParams,referenceLMS,comparisonLMS)
%
% Syntax:
%     opponentContrast = LMSToOpponentContrast(colorDiffParams,referenceLMS,comparisonLMS)
%
% Description:
%     Convert LMS representation to opponent contrast respresentation
%
%     When the type field is 'opponentColor', converts comparison to
%     contrast with respect to reference, then transforms to a lum, rg, by
%     opponent representation. The contrast differences in each direction
%     are scaled by the weights specified in the parameters structure.
%
% Inputs:
%     colorDiffParams       - Structure of parameters describing color
%                             difference model, with fields:
%                               type: String with model type.  Only current
%                                 option is 'opponentContrast', but you never
%                                 know what the future will bring.
%                               LMRatio: LM cone ratio, determines opponent
%                                 luminance sensitivity.
%                               lumWeight: Weight on luminance mechanism.
%                               rgWeight: Weight on rg mechanism.
%                               byWeight: Weight on by mechanism.
%     referenceLMS           - LMS coordinates of the reference with
%                              respect to which contrast is computed.
%     comparisonLMS          - LMS coordinates of the stimulus whose
%                              contrast is computed.
%
% Outputs:
%   opponentContrast         - Contrast representation of the comparison
%
% Optional key/value pairs:
%   None:
%
% See also: OpponentContrastToLMS, ComputeMatchDiff
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
    opponentContrast = LMSToOpponentContrast(colorDiffParams,referenceLMS,2*referenceLMS)
    if (any(abs(opponentContrast - [1 0 0]')) > 1e-6)
       error('Don''t get right answer in simple test case.');
    end
%}

switch (colorDiffParams.type)
    case 'opponentContrast'
        
        % Build up matrix
        M = diag([colorDiffParams.lumWeight colorDiffParams.rgWeight colorDiffParams.byWeight])*[ [colorDiffParams.LMRatio 1 0]/(colorDiffParams.LMRatio+1) ; ...
              [1 -1 0] ; ...
              [-0.5 -0.5 1] ];
        
        % Go to cone contrast         
        coneContrast = (comparisonLMS - referenceLMS) ./ referenceLMS;
        
        % And then to opponent contrast
        opponentContrast = M*coneContrast;
        
    otherwise
        error('Unknown type field passed in colorDiffParams structure');
end

