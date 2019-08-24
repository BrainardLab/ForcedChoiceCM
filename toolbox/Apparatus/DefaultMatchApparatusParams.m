function matchApparatusParams = DefaultMatchApparatusParams(type,S)
% Generate default color difference params structure
%
% Syntax:
%    matchApparatusParams = DefaultApparatusParams(type,S)
%
% Description:
%    Generate a structure describing the properties of the apparatus, and
%    also that holds information about the stimulus parameters.
%
%    The input string type allows some flexibility about the description.
%
% Inputs:
%     type                          - String specifying cone parameterization type.
%                                     'monochromatic': Primary lights are monochromatic.
%     S                             - Wavelength sampling information to use
%
% Outputs:
%     matchApparatusParams          - Structure with field for each parameter.
%
% Optional key/value pairs:
%     None.
%
% See also: StimulusVecToParams, StimulusParamsToVec
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    S = [400 1 301];
    apparatusParams = DefaultMatchApparatusParams('monochromatic',S);
    apparatusParams
%}

%% Set type
matchApparatusParams.type = type;

%% Do the right thing according to type
switch (matchApparatusParams.type)
    case 'monochromatic'
                
        % Wavelength sampling
        matchApparatusParams.S = S;
        matchApparatusParams.wls = SToWls(matchApparatusParams.S);
        
        % Apparatus parameters
        matchApparatusParams.primaryWavelength1 = 430;
        matchApparatusParams.primaryWavelength2 = 545;
        matchApparatusParams.primaryWavelength3 = 670;
        
        % Compute indices and set spectra below
        matchApparatusParams.primaryIndex1 = find(matchApparatusParams.wls == matchApparatusParams.primaryWavelength1);
        matchApparatusParams.primaryIndex2 = find(matchApparatusParams.wls == matchApparatusParams.primaryWavelength2);
        matchApparatusParams.primaryIndex3 = find(matchApparatusParams.wls == matchApparatusParams.primaryWavelength3);
        
        matchApparatusParams.unitPrimarySpectrum1 = zeros(size(matchApparatusParams.wls));
        matchApparatusParams.unitPrimarySpectrum1(matchApparatusParams.primaryIndex1) = 1;
        matchApparatusParams.unitPrimarySpectrum2 = zeros(size(matchApparatusParams.wls));
        matchApparatusParams.unitPrimarySpectrum2(matchApparatusParams.primaryIndex2) = 1;
        matchApparatusParams.unitPrimarySpectrum3 = zeros(size(matchApparatusParams.wls));
        matchApparatusParams.unitPrimarySpectrum3(matchApparatusParams.primaryIndex3) = 1;
        
        matchApparatusParams.primaryBasis = [matchApparatusParams.unitPrimarySpectrum1 ... 
            matchApparatusParams.unitPrimarySpectrum2 ... 
            matchApparatusParams.unitPrimarySpectrum3];
        
    otherwise
        error('Unknown apparatus parameters type passed.');
end

