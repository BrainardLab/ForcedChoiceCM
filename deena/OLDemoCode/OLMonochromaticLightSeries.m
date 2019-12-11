function OLMonochromaticLightSeries(varargin)
% Presents a series of monochromatic lights on OneLight box
%
% Syntax:
%t   OLMonochromaticLightSeries
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
%
% Optional key-value pairs:
%    'checkPlots'   - logical indicating whether to plot predicted and
%                     target spds (mostly for initial error checking).
%                     Default is false

% History:
%    9/10/19   dce    Wrote program

% Parse input
p = inputParser;
p.addParameter('checkPlots', false, @(x) (islogical(x)));
p.parse(varargin{:});

% Get the calibration structure
cal = OLGetCalibrationStructure;

% Define target wavelengths and key parameters. Our target wavelengths
% range from 490-690nm in 10nm increments.
wls = 490:10:690;
fullWidthHalfMax = 20;
lambda = 0.001;
delaySecs = 0.3; % time in seconds that each wavelength is displayed
wlsPos = 1; % Initial position in wavelength array

% Initialize arrays for storing data
numWls = length(wls);
numCols = cal.describe.numColMirrors; % number of mirrors in each OL column
[spdLength,primaryLength,~] =  size(cal.computed.pr650M);

targetSpds = zeros(spdLength,numWls);
primaries = zeros(primaryLength,numWls);
settings = zeros(primaryLength,numWls);
startStops = zeros(numWls,2,numCols); % one column each for start and stop

% Convert wavelengths to OneLight settings 
for i = 1:numWls
    targetSpds(:,i) = OLMakeMonochromaticSpd(cal, wls(i), fullWidthHalfMax)/3; % where does 3 come from?
    primaries(:,i) = OLSpdToPrimary(cal, targetSpds(:,i), 'lambda', lambda);
    settings(:,i) = OLPrimaryToSettings(cal, primaries(:,i)); % gamma correction
    [start,stop] = OLSettingsToStartsStops(cal, settings(:,i));
    startStops(i,1,:) = start;
    startStops(i,2,:) = stop;
end

% Plot target and predicted spds to check that the OneLight can produce the
% values we want. This option is turned off by default.
if p.Results.checkPlots
    predictedSpds = zeros(201,numWls);
    for i = 1:numWls
        predictedSpds(:,i) = OLPrimaryToSpd(cal,primaries(:,i));
        figure(i); hold on
        plot(SToWls(cal.describe.S),targetSpds(:,i),'r');
        plot(SToWls(cal.describe.S),predictedSpds(:,i),'b');
        title(['target and predicted spds: ', num2str(wls(i)), ' nm'])
        xlabel('wavelength')
        ylabel('intensity')
        legend('Target spd','Predicted spd')
        hold off
    end
end

% Enable character listening
mglDisplayCursor(0);
ListenChar(2);
FlushEvents;
fprintf('Starting display loop \n'); 

% Open OneLight
ol = OneLight;

% Loop through lights until the user presses a key
while(~CharAvail)
    % Display chosen wavelength on OneLight
    ol.setMirrors(squeeze(startStops(wlsPos,1,:))', squeeze(startStops(wlsPos,2,:))');
    
    fprintf('Waiting %g seconds for wavelength %g \n',delaySecs,wls(wlsPos));
    mglWaitSecs(delaySecs); % time delay
    
    % increment position in wavelength ar
    wlsPos = wlsPos + 1;
    if wlsPos > numWls
        wlsPos = 1;
    end
end

% Clean up once user exits
fprintf('User exited the program \n'); 
ListenChar(0);
mglDisplayCursor(1);
end