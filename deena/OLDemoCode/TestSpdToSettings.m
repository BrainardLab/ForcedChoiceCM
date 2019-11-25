% Script to try to help debug why OLSpdToSettings does odd things at edge
% of spectral range

% Get calibration information
cal = OLGetCalibrationStructure;

% Calculate sizes and initialize spd arrays
wls = cal.computed.pr650Wls;
S = WlsToS(wls);
fullWidthHalfMax = 20;
lambda = 0.001;
testWl = 580;

% Initialize arrays
testSpdNominal = zeros(length(wls), 1);
testSpdPredicted = zeros(length(wls), 1);

primaryScale = 1;
testSpdNominal = OLMakeMonochromaticSpd(cal, testWl, fullWidthHalfMax)*primaryScale;
[testPrimary,testSpdPredicted] = OLSpdToPrimary(cal, testSpdNominal, 'lambda', lambda, ...
    'whichSpdToPrimaryMin', 'leastSquares', ...
    'verbose', true);

figure;
hold on
plot(wls,testSpdPredicted,'r','LineWidth',4);
plot(wls,testSpdNominal,'k','LineWidth',2);
