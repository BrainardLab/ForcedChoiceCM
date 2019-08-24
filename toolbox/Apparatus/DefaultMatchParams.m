function matchParams = DefaultMatchParams(matchApparatusParams)
% Generate default matching light structure
%
% Syntax:
%    matchParams = DefaultMatchParams(matchApparatusParams)
%
% Description:
%    Generate a structure describing the properties of the matching light.
%
%    The input string type allows passing information about the appartus.
%
% Inputs:
%     matchApparatusParams          - Structure with field for each parameter.
%
% Outputs:
%     matchParams                   - Structure with field for each parameter.
%
% Optional key/value pairs:
%     None.
%
% See also: DefaultMatchApparatusParams, SetMatchParams, DefaultTestParams, SetTestParams, StimulusVecToParams, StimulusParamsToVec
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    S = [400 1 301];
    matchApparatusParams = DefaultMatchApparatusParams('monochromatic',S);
    matchParams = DefaultMatchParams(matchApparatusParams);
    matchParams
%}

%% Do the right thing according to type
switch (matchApparatusParams.type)
    case 'monochromatic'
        
        % Apparatus parameters
        matchParams.type = matchApparatusParams.type;
        matchParams.primary = [0.5 0.5 0.5]';
        matchParams.spectrum = matchApparatusParams.primaryBasis*matchParams.primary;
        
    otherwise
        error('Unknown apparatus parameters type passed.');
end

