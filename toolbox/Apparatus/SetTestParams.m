function testParams = SetTestParams(testParams)
% Make test parameter structure match current setable values
%
% Syntax:
%    testParams = SetTestParams(testParams)
%
% Description:
%    After the user mucks with certain fields, the structure needs to be
%    reinitialized. This routine does it.
%
%    If I were a better person I'd make this an object and have all this
%    happen in the slick more bullet proof object oriented way.
%
% Inputs:
%     testParams                     - Structure with field for each parameter.
%
% Outputs:
%     testParams                     - Structure with field for each parameter.
%
% Optional key/value pairs:
%     None.
%
% See also: StimulusVecToParams, StimulusParamsToVec, DefaultTestParams, DefaultMatchApparatusParams
%

% History:
%   08/10/19  dhb  Wrote it.

% Examples:
%{
    S = [400 1 301];
    testParams = DefaultTestParams('monochromatic',S)
    testParams.testWavelength = 630;
    testParams = SetTestParams(testParams)
%}

%% Dp the right thing according to type
switch (testParams.type)
    case {'monochromatic', 'rayleigh'}
                
        % Compute indices and set spectra below
        testParams.testIndex = find(testParams.wls == testParams.testWavelength);
        testParams.unitTestSpectrum = zeros(size(testParams.wls));
        testParams.unitTestSpectrum(testParams.testIndex) = 1;
        
    otherwise
        error('Unknown apparatus parameters type passed.');
end

