function OLMonochromaticLightSeries
% Presents a series of monochromatic lights on OneLight box
%
% Syntax:
%   OLMonochromaticLightSeries
%
% Description:
%    Sends commands to the OneLight box to present a series of
%    monochromatic lights. The lights presented range from 490-690nm in
%    10nm increments. Listens for user input and exits when a key is
%    pressed.
%
% Inputs
%    None
%
% Outputs:
%    None

% History:
%    9/10/19   dce    Wrote program

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Define target wavelengths and key parameters. Our target wavelengths
% range from 490-690nm in 10nm increments.
wls = 490:10:690;
fullWidthHalfMax = 20;
lambda = 0.001;
wlsPos = 1; % Initial position in wavelength array 

% Add conversion from wls to settings

% Enable character listening
mglDisplayCursor(0);
ListenChar(2);
FlushEvents;

% Circle through lights until the user presses a key 
while(~CharAvail)
    % Calculate start and stop points and display on OneLight
    [starts1, stops1] = OLSettingsToStartsStops(cal, settings1);
    ol = OneLight;
    ol.setMirrors(starts1, stops1);
    % add delay
    % increment position for next wavelength to show
    wlsPos = wlsPos + 1;
    if wlsPos == length(wls)
        wlsPos = 1;
    end
end

% Clean up
ListenChar(0);
mglDisplayCursor(1);
end