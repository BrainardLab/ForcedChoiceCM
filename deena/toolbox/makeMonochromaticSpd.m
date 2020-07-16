function spd = makeMonochromaticSpd(wavelength,scaleFactor,S)
% Syntax:
%   makeMonochromaticSpd(wavelength,scaleFactor,S)
%
% Description:
%    Constructs a monochromatic spd that has a value given by scaleFactor
%    at the peak wavelength, and 0 elsewhere.

% Inputs:
%    wavelength  -Integer wavelength for spd peak.
%    scaleFactor -Number between 0 and 1 for scaling the spd height.
%    S           -Wavelength sampling for generating the spd, in the form 
%                 [start delta nTerms].
%
% Outputs:
%    spd         -Monochromatic spd 
%
% Optional key-value pairs:
%    None

% History
%    dce    7/16/20   -Moved from computePredictedRayleighMatch

wls = SToWls(S);
spd = zeros(size(wls));
spd(wls==wavelength) = scaleFactor;
end