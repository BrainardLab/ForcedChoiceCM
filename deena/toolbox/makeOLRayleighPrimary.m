function spd = makeOLRayleighPrimary(wavelength,varargin)
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
%    wavelength     -Integer wavelength for spd peak.
%
% Outputs:
%    spd            -Column vector of the narrowband spd.
%
% Optional key-value pairs:
%    scaleFactor   -Number between 0 and 1 that scales the height of the
%                   spd. Default is 1.
%    nominal       -Logical indicating to return nominal spds rather than
%                   predicted values. Default is false.
%    subtractDark  -Logical indicating to subtract the dark spd from the
%                   final predicted spd (if it is not nominal). Defult is 
%                   true.

% History
%    dce    07/08/20  - Wrote it
%    dce    07/14/20  - Got rid of division by 3 for p2
%    dce    07/16/20  - Added scaling to this file
%    dce    07/17/20  - Made subtractDark optional

% Parse input
p = inputParser;
p.addParameter('scaleFactor',1,@(x)(isnumeric(x)));
p.addParameter('nominal',false,@(x)(islogical(x)));
p.addParameter('subtractDark',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Base parameters
fullWidthHalfMax = 20;
lambda = 0.001;

% Load default calibration
cal = OLGetCalibrationStructure('CalibrationType',...
    getpref('ForcedChoiceCM','currentCal'));
[~,settingsLength] = size(cal.computed.pr650M);
darkSpd = OLPrimaryToSpd(cal,zeros(settingsLength,1));

if wavelength == 0
    spd = darkSpd;
else
    % Find spd
    spdRel = OLMakeMonochromaticSpd(cal,wavelength,fullWidthHalfMax);
    spdNominal = OLMaximizeSpd(cal,spdRel,'lambda',lambda)...
        *p.Results.scaleFactor;
    
    % Convert to predicted spd
    if p.Results.nominal
        spd = spdNominal;
    else
        [~,~,spd] = OLSpdToSettings(cal,spdNominal+darkSpd,'lambda',...
            lambda);
        if p.Results.subtractDark
            spd = spd-darkSpd;
        end 
    end
end

% Save results
if wavelength==0
    file = sprintf('OLRayleighPrimary_dark.mat');
elseif p.Results.nominal 
    file = sprintf('OLRayleighPrimary_%g_%g_nominal.mat',wavelength,...
        p.Results.scaleFactor);
else 
    file = sprintf('OLRayleighPrimary_%g_%g_predicted.mat',wavelength,...
        p.Results.scaleFactor);
end 
save(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedSpds', file),'spd');
end