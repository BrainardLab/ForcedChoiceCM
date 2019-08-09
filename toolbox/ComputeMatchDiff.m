function colorDiff = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%
% Syntax:
%     colorDiff = ComputeMatchDiff(colorDiffParams,referenceLMS,comparisonLMS)
%
% Description:
%     Compute a single number color difference between two vectors
%     specified in cone coordinates.
%
% Inputs:
%
% Outputs:
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

