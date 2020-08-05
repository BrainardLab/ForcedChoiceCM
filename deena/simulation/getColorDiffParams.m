function colorDiffParams = getColorDiffParams(opponentParams)
% Generate color difference params structure with weights specified by
% input
%
% Syntax:
%    colorDiffParams = getColorDiffParams(type)
%
% Description:
%    Generate a structure describing and observer's cone fundamentals, with
%    reasonable defaults.
%
%    The input string type allows some flexibility about the description.
%
% Inputs:
%     opponentParams    - 4-element numeric vector with opponent contrast
%                         parameters. (1) specifies luminance weight, (2) 
%                         specifies rg weight, (3) specifies by weight, and 
%                         (4) specifies noise sd.
% Outputs:
%     colorDiffParams   - Structure with field for each parameter.
%
% Optional key/value pairs:
%     None.
%
% See also: ObserverVecToParams, ObserverParamsToVec
%

% History:
%   08/05/20  dce  Adapted from DefaultConeParams.

colorDiffParams = struct();
colorDiffParams.type = 'opponentContrast';
colorDiffParams.LMRatio = 2;
colorDiffParams.lumWeight = opponentParams(1);
colorDiffParams.rgWeight = opponentParams(2);
colorDiffParams.byWeight = opponentParams(3);
colorDiffParams.noiseSd = opponentParams(4);
colorDiffParams.M = GetOpponentContrastMatrix(colorDiffParams);
end

