function apparatusParams = DefaultApparatusParams(type)
% Generate default color difference params structure
%
% Syntax:
%    apparatusParams = DefaultApparatusParams(type)
%
% Description:
%    Generate a structure describing the properties of the apparatus, and
%    also that holds information about the stimulus parameters.
%
%    The input string type allows some flexibility about the description.
%
% Inputs:
%     type                          - String specifying cone parameterization type.
%                                     'monochromatic': Test and primary
%                                      lights are monochromatic.
% Outputs:
%     apparatusParams               - Structure with field for each parameter.
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
    apparatusParam = DefaultApparatusParams('monochromatic');
    apparatusParam
%}

%% Set type
apparatusParams.type = type;

%% Dp the right thing according to type
switch (apparatusParams.type)
    case 'monochromatic'
                
        % Wavelength sampling
        apparatusParams.S = [400 1 301];
        apparatusParams.wls = SToWls(apparatusParams.S);
        
        % Apparatus parameters
        apparatusParams.testWavelength = 520;
        apparatusParams.primaryWavelength1 = 430;
        apparatusParams.primaryWavelength2 = 545;
        apparatusParams.primaryWavelength3 = 670;
        
        % Compute indices and set spectra below
        apparatusParams.testIndex = find(apparatusParams.wls == apparatusParams.testWavelength);
        apparatusParams.primaryIndex1 = find(apparatusParams.wls == apparatusParams.primaryWavelength1);
        apparatusParams.primaryIndex2 = find(apparatusParams.wls == apparatusParams.primaryWavelength2);
        apparatusParams.primaryIndex3 = find(apparatusParams.wls == apparatusParams.primaryWavelength3);
        apparatusParams.unitTestSpectrum = zeros(size(apparatusParams.wls));
        apparatusParams.unitTestSpectrum(apparatusParams.testIndex) = 1;
        apparatusParams.unitPrimarySpectrum1 = zeros(size(apparatusParams.wls));
        apparatusParams.unitPrimarySpectrum1(apparatusParams.primaryIndex1) = 1;
        apparatusParams.unitPrimarySpectrum2 = zeros(size(apparatusParams.wls));
        apparatusParams.unitPrimarySpectrum2(apparatusParams.primaryIndex2) = 1;
        apparatusParams.unitPrimarySpectrum3 = zeros(size(apparatusParams.wls));
        apparatusParams.unitPrimarySpectrum3(apparatusParams.primaryIndex3) = 1;
        apparatusParams.primaryBasis = [apparatusParams.unitPrimarySpectrum1 ... 
            apparatusParams.unitPrimarySpectrum2 ... 
            apparatusParams.unitPrimarySpectrum3];
        
    otherwise
        error('Unknown apparatus parameters type passed.');
end

