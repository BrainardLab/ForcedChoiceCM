function spd = makeOLRayleighPrimary(wavelength,varargin)
% Generates a narrowband spd for use in OneLight Rayleigh calculations
%
% Syntax:
%   makeOLRayleighPrimary(wavelength)
%
% Description:
%    Takes in a peak wavelength as input, then constructs an associated spd
%    based on the default OneLight calibration. If wavelength of 0 is 
%    entered, the program constructs a dark spd. The spd
%    is saved in a file with the naming convention
%    'OLRayleighPrimary_wavelength_nominal/predicted.mat', where it can be 
%    accessed for future use.
%
%    When predicted spds are used, two additional options exist as key-
%    value pairs. 'subtractDark' subtracts the dark spd before returning
%    the final spd. Additionally, 'subtractAroundWl' replaces a portion of
%    the predicted spd around the specified wavelength with the nominal
%    spd. Both of these options should be used when producing spectra to
%    use with computePredictedRayleighMatch. These options are reflected in
%    the file-naming convention when used. Note that the options are
%    ignored if 'nominal' is used.

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
%    subtractAroundWl  -Interger wavelength to subtract around. When a
%                       nonzero value is specified, replaces the predicted 
%                       values of the spd around that range with the 
%                       nominal values. Default is 0.

% History
%    dce    07/08/20  - Wrote it
%    dce    07/14/20  - Got rid of division by 3 for p2
%    dce    07/16/20  - Added scaling to this file
%    dce    07/17/20  - Made subtractDark optional
%    dce    07/21/20  - Made option to subtract around a wavelength, and 
%                       changed file-naming convention 

% Parse input
p = inputParser;
p.addParameter('scaleFactor',1,@(x)(isnumeric(x)));
p.addParameter('nominal',false,@(x)(islogical(x)));
p.addParameter('subtractDark',true,@(x)(islogical(x)));
p.addParameter('subtractAroundWl',0,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Base parameters
fullWidthHalfMax = 20;
lambda = 0.001;

% Load default calibration
cal = OLGetCalibrationStructure('CalibrationType',...
    getpref('ForcedChoiceCM','currentCal'));
wls = cal.computed.pr650Wls;
[~,settingsLength] = size(cal.computed.pr650M);
darkSpd = OLPrimaryToSpd(cal,zeros(settingsLength,1));

if wavelength == 0
    spd = darkSpd;
else
    % Find spd
    spdRel = OLMakeMonochromaticSpd(cal,wavelength,fullWidthHalfMax);
    spdNominal = OLMaximizeSpd(cal,spdRel,'lambda',lambda)...
        *p.Results.scaleFactor;
    
    % Convert to predicted spd, if desired
    if p.Results.nominal
        spd = spdNominal;
    else  
        [~,~,spd] = OLSpdToSettings(cal,spdNominal+darkSpd,'lambda',...
            lambda);
        % If you're subtracting around a wavelength, define the limit 
        % indices and replace the predicted spectrum with the nominal
        % spectrum in that portion
        if p.Results.subtractAroundWl~=0
           wl1 = p.Results.subtractAroundWl-fullWidthHalfMax;
           if wl1 < 0 
               wl1 = 0;
           end 
           wl2 = p.Results.subtractAroundWl+fullWidthHalfMax;
           if wl2 > wls(end)
               wl2 = wls(end);
           end 
           ind1 = find(wls==wl1);          
           ind2 = find(wls==wl2);
           spd(ind1:ind2) = spdNominal(ind1:ind2)+darkSpd(ind1:ind2);
        end 
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
    baseFile = sprintf('OLRayleighPrimary_%g_%g_predicted',wavelength,...
        p.Results.scaleFactor);
    if p.Results.subtractDark
        baseFile = [baseFile '_subtractDark'];
    end
    if p.Results.subtractAroundWl ~= 0 
        wlName = sprintf('_subtract%g',p.Results.subtractAroundWl);
        baseFile = [baseFile wlName];
    end 
    file = [baseFile '.mat'];
end 
save(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedSpds', file),'spd');
end