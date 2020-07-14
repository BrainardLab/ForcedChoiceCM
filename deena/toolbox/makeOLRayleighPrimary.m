function spd = makeOLRayleighPrimary(wavelength)
% Generates a narrowband spd for use in OneLight Rayleigh calculations
%
% Syntax:
%   makeOLRayleighPrimary(wavelength)
%
% Description:
%    Takes in a peak wavelength as input, then constructs an associated spd
%    based on the default OneLight calibration. If 'p2' is selected, the
%    spd is divided by 3 before being scaled to its maximum. If a
%    wavelength of 0 is entered, the program constructs a dark spd. The spd
%    is saved in a file with the naming convention
%    'OLRayleighPrimary_wavelength.mat', where it can be accessed for
%    future use.

% Inputs:
%    wavelength  -Integer wavelength for spd peak.
%
% Outputs:
%    none (saves a file)
%
% Optional key-value pairs:

% History
%    dce    07/08/20  - Wrote it
%    dce    07/14/20  - Got rid of division by 3 for p2

% Load default calibration
cal = OLGetCalibrationStructure('CalibrationType',getpref('ForcedChoiceCM','currentCal'));

% Make spd
if wavelength == 0 % Make dark spd
    [~,settingsLength] = size(cal.computed.pr650M);
    spd = OLPrimaryToSpd(cal,zeros(settingsLength,1));
    
else               % Make narrowband spd
    % Base parameters
    fullWidthHalfMax = 20;
    lambda = 0.001;
    
    % Calculation
    spdRel = OLMakeMonochromaticSpd(cal,wavelength,fullWidthHalfMax);
    spd = OLMaximizeSpd(cal,spdRel,'lambda',lambda);
end

% Save results
file = sprintf('OLRayleighPrimary_%g.mat',wavelength);
save(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedSpds', file),'spd');
end