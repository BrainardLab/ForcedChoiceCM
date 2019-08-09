function colorDiff = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%
% Syntax:
%     colorDiff = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%
% Description:
%     Compute a single number color difference between two vectors
%     specified in cone coordinates. 
%
%     When the type field is 'opponentColor', converts comparison to
%     contrast with respect to reference, then transforms to a lum, rg, by
%     opponent representation and computes vector length.  The differences
%     in each opponent mechanism are scaled by the weight parameters before
%     vector length is taken.
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
%
% Outputs:
%   colorDiff                - Single number color difference measure.
%
% Optional key/value pairs:
%   None:
%
% See also:
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

switch (colorDiffParams.type)
    case 'opponentContrast'
        
        % Convert to LMS contrast
        coneContrast = (comparisonLMS-referenceLMS)./referenceLMS;
        
        % Approximate three post-receptoral constrasts
        lumContrast = (colorDiffParams.LMRatio*coneContrast(1)+coneContrast(2))/(colorDiffParams.LMRatio+1);
        rgContrast = coneContrast(1)-coneContrast(2);
        byContrast = coneContrast(3) - (coneContrast(1)+coneContrast(2))/2;
        
        % Take appropriate weighted vector length
        colorDiff = sqrt((colorDiffParams.lumWeight*lumContrast)^2 + (colorDiffParams.rgWeight*rgContrast)^2 + (colorDiffParams.byWeight*byContrast)^2);
        
    otherwise
        error('Unknown type field passed in colorDiffParams structure');
end

