function matchParams = SetMatchParams(matchApparatusParams,matchParams)
% Make match parameter structure match current setable values
%
% Syntax:
%    matchParams = SetMatchParams(matchApparatusParams,matchParams)
%
% Description:
%    After the user mucks with certain fields, the structure needs to be
%    reinitialized. This routine does it.
%
%    If I were a better person I'd make this an object and have all this
%    happen in the slick more bullet proof object oriented way.
%
% Inputs:
%     matchApparatusParams           - Structure describing match apparatus.
%     matchParams                     - Structure describing match stimulus.
%
% Outputs:
%     matchParams                     - Structure describing match stimulus.
%
% Optional key/value pairs:
%     None.
%
% See also: DefaultMatchParams, StimulusVecToParams, StimulusParamsToVec, DefaultTestParams, DefaultMatchApparatusParams
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    S = [400 1 301];
    matchApparatusParams = DefaultMatchApparatusParams('monochromatic',S);
    matchParams = DefaultMatchParams(matchApparatusParams)
    matchParams.primary = [0.2 0.7 0.3]';
    matchParams = SetMatchParams(matchApparatusParams,matchParams)
%}

% Check
if (~strcmp(matchParams.type,matchApparatusParams.type))
    error('Type mismatch between apparatus and stimulus structures');
end

%% Do the right thing according to type
switch (matchParams.type)
    case 'monochromatic'
        % Update
        matchParams.spectrum = matchApparatusParams.primaryBasis*matchParams.primary;

    otherwise
        error('Unknown apparatus parameters type passed.');
end

