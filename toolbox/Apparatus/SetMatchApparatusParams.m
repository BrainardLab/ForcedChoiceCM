function matchApparatusParams = SetMatchApparatusParams(matchApparatusParams)
% Make match apparaus parameter structure match current setable values
%
% Syntax:
%    matchApparatusParams = SetMatchApparatusParams(matchApparatusParams)
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
%
% Outputs:
%     matchApparatusParams           - Updated structure describing match stimulus.
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
    matchApparatusParams = DefaultMatchApparatusParams('rayleigh',S);
    matchApparatusParams.primaryWavelength1 = 490;
    matchApparatusParams = SetMatchApparatusParams(matchApparatusParams)
%}

%% Do the right thing according to type
switch (matchApparatusParams.type)
    case 'monochromatic'
        error('Have not implemented this for type ''monochrmatic'' yet');

    case 'rayleigh'
                  
        % Compute indices and set spectra below
        matchApparatusParams.primaryIndex1 = find(matchApparatusParams.wls == matchApparatusParams.primaryWavelength1);
        matchApparatusParams.primaryIndex2 = find(matchApparatusParams.wls == matchApparatusParams.primaryWavelength2);
        
        matchApparatusParams.unitPrimarySpectrum1 = zeros(size(matchApparatusParams.wls));
        matchApparatusParams.unitPrimarySpectrum1(matchApparatusParams.primaryIndex1) = 1;
        matchApparatusParams.unitPrimarySpectrum2 = zeros(size(matchApparatusParams.wls));
        matchApparatusParams.unitPrimarySpectrum2(matchApparatusParams.primaryIndex2) = 1;
        
        matchApparatusParams.primaryBasis = [matchApparatusParams.unitPrimarySpectrum1 ... 
            matchApparatusParams.unitPrimarySpectrum2];

    otherwise
        error('Unknown apparatus parameters type passed.');
end

