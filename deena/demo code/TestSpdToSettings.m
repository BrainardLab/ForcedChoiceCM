% Show how to automatically find largest scale factor that
% works reasonably with any desired test or primary.

% History:
%   11/27/19  dhb  First more or less working version
%   12/24/19  dhb  Edit to use new routine OLMaximizeSpd

%% Clear
clear; close all;

%% Get OL calibration information
cal = OLGetCalibrationStructure;

%% Parameters
fullWidthHalfMax = 20;
lambda = 0.001;
testWl = 500;
errorTolerance = 0.05;
maxScaleFactor = 10;
nScaleFactors = 100;

%% Make nominal spd
wls = cal.computed.pr650Wls;
S = WlsToS(wls);
testSpdIncrRel = OLMakeMonochromaticSpd(cal, testWl, fullWidthHalfMax);

%% Get the "dark" spd, that which comes out when primarys are at 0
darkSpd = OLPrimaryToSpd(cal,zeros(size(cal.computed.D,2),1));

%% Maximize the spd
[maxTestSpdIncr,spdScaleFactor] = OLMaximizeSpd(cal,testSpdIncrRel, ...
    'lambda', lambda, 'showPlot', true);

%% Find primaries for a number scale factors
%
% Compute relative error for each.  Relative error
% is used because otherwise it will increase with
% scale factor even if the shape is the same.
%
% Differential mode looks at the effect of the primaries
% relative to the light that comes out when the input
% is set to 0. For our purposes, this is what we want.
primaryScales = linspace(0,1,nScaleFactors);
for ii = 1:length(primaryScales)
    primaryScale = primaryScales(ii);
    [testPrimary,maxTestSpdIncrPredicted] = OLSpdToPrimary(cal, primaryScale*maxTestSpdIncr+darkSpd, 'lambda', lambda, ...
        'whichSpdToPrimaryMin', 'leastSquares', ...
        'verbose', false);
end

%% Use the max scale factor and plot nominal and predicted spds
[testPrimary,maxTestSpdIncrPredicted] = OLSpdToPrimary(cal, maxTestSpdIncr+darkSpd, 'lambda', lambda, ...
        'whichSpdToPrimaryMin', 'leastSquares', ...
        'verbose', false);
    
% The dark spd is added back in here, so that we represent what
% we actually think will reach the eye.
figure;
hold on
plot(wls,maxTestSpdIncrPredicted,'r','LineWidth',4);
plot(wls,maxTestSpdIncr+darkSpd,'k','LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Power');
legend({'Predicted', 'Nominal'});

%% Show how to get settings and starts/stops
primarySettings = OLPrimaryToSettings(cal,testPrimary);
