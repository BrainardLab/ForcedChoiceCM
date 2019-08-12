function testParams = DefaultTestParams(type,S)
% Generate default test light params structure
%
% Syntax:
%    testParams = DefaultTestParams(type,S)
%
% Description:
%    Generate a structure describing the properties of the test light.
%
%    The input string type allows some flexibility about the description.
%
% Inputs:
%     type                          - String specifying cone parameterization type.
%                                     'monochromatic': Test lights are monochromatic.
%     S                             - Wavelength sampling information to use
%
% Outputs:
%     testParams                    - Structure with field for each parameter.
%
% Optional key/value pairs:
%     None.
%
% See also: SetTestParams, StimulusVecToParams, StimulusParamsToVec, DefaultMatchApparatusParams
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    S = [400 1 301];
    testParams = DefaultTestParams('monochromatic',S);
    testParams
%}

%% Set type
testParams.type = type;

%% Do the right thing according to type
switch (testParams.type)
    case 'monochromatic'
                
        % Wavelength sampling
        testParams.S = S;
        testParams.wls = SToWls(testParams.S);
        
        % Apparatus parameters
        testParams.testWavelength = 520;
        
        % Compute indices and set spectra below
        testParams.testIndex = find(testParams.wls == testParams.testWavelength);
        testParams.unitTestSpectrum = zeros(size(testParams.wls));
        testParams.unitTestSpectrum(testParams.testIndex) = 1;
        testParams.testIntensity = 1;
        
    otherwise
        error('Unknown apparatus parameters type passed.');
end

