function [testAdjustedSpd,primaryMixtureSpd,testIntensity,lambda] = ...
    computePredictedRayleighMatch(p1Spd,p2Spd,testSpd,observer,varargin)
% Calculates a predicted Rayleigh match analytically
%
% Syntax:
%    [testAdjustedSpd,primaryMixtureSpd,testIntensity,lambda] = ...
%        computePredictedRayleighMatch(p1Spd,p2Spd,testSpd,observerParams, ...
%        opponentParamsVec)
%
% Description:
%    Computes a Rayleigh match analytically when given the associated
%    primary and test wavelengths and the individual difference parameters
%    for the observer (based on the Asano model). Makes a monochromatic spd
%    for each of the wavelengths, then solves a matrix equation to find the
%    primary ratio and test intensity which would lead to a match. If a
%    noise standard deviation is entered for the observer, then Gaussian
%    noise is added to the lambda and testIntensity parameters. There is
%    also an option to plot the chosen match spds and their cone effects.

%    Unlike searchPredictedRayleighMatch, which finds the best available
%    match, this program assumes that the primary spectra are perfect
%    linear combinations of p1 and p2. This assumption holds for nominal
%    OneLight spectra, but not predicted. The linearity assumption also
%    breaks down somewhat when adjustmentLength is small. 
%
% Inputs:
%    p1Spd           -nx1 spd of the first primary light.
%    p2Spd           -nx1 of the second primary light.
%    testSpd         -nx1 of the test light.
%    observer        -Structure of a simulated observer (see
%                     genRayleighObserver for details)
%
% Outputs:
%    primaryMixtureSpd    -Predicted primary spd for the match
%    testAdjustedSpd      -Predicted test spd for the match
%    lambda               -Scale factor for primary 1 proportion, between 0
%                          and 1.
%    testIntensity        -Scale factor for test intensity, between 0 and 1
%
% Optional key-value pairs:
%    'addDarkSpd'      -Logical.  Default is true. When true, adds a dark
%                       spd as specified by 'darkSpd'.
%    'darkSpd'         -nx1 dark spd vector, default is []
%                       In case of [] and monochromatic is false, the dark
%                       light from current OL calibration is added to
%                       spectra before returning. Send in an nx1 vector of
%                       zeros to prevent this behavior.
%    'S'               -Wavelength sampling for cone calculations, in the
%                       form [start delta nTerms]. Default is [380 2 201],
%                       which is the OneLight convention. Note that S can
%                       only vary from the OneLight convention when
%                       monochromatic lights are used. Otherwise, the
%                       routine will throw an error.

% History
%    dce    6/29/20   -Adapted from example code from dhb
%    dce    7/01/20   -Added noise
%    dce    7/03/20   -Add noise to lambda and testIntensity, not cone res
%    dce    7/06/20   -Added option to use OL spds, added scale factors,
%                      made no noise the default
%    dce    7/14/20   -Added predicted matching option
%    dce    7/16/20   -Separated spd creation into a separate file
%    dce    7/22/20   -Added option to pass in dark spd
%    dhb    08/08/20  -Remame opponentParams -> opponentParamsVec for clarity.
%    dce    08/10/20  -Rename monochromatic key-value pair, changed to take
%                      in observer struct instead of params vector.

%% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('addDarkSpd',true,@(x)(islogical(x)));
p.addParameter('darkSpd',[],@(x)(isnumeric(x)));
p.parse(varargin{:});

%% Find observer cone responses to the three base spds. 
% Each vector contains two elements-one for the L cone, one for the M cone.
T_LM = observer.T_cones(1:2,:);         % L and M cone fundamentals
p1Res = T_LM*p1Spd;
p2Res = T_LM*p2Spd;
testRes = T_LM*testSpd;

%% Set up a matrix equation to solve for the optimal primary ratio and test intensity
%
% We want this to be true
%   T_LM*(testIntensity)*testSpd == T_LM*(lambda*primary1Spd + (1-lambda)*primary2Spd)
%   0 = T_LM*[(lambda*primary1Spd + (1-lambda)*primary2Spd) - (testIntensity)*testSpd)
%   0 = T_LM*primary1Spd * lambda + T_LM*primary2Spd * (1-lambda) - T_LM*testSpd * testIntensity)

% So
%    0 = p1Res(1) * lambda + p2Res(1) - p2Res(1)*lambda - testRes(1)*testIntensity
%    0 = p1Res(2) * lambda + p2Res(2) - p2Res(2)*lambda - testRes(2)*testIntensity
%
%    -p2Res(1) = [p1Res(1)- p2Res(1)]*lambda - testRes(1)*testIntensity
%    -p2Res(2) = [p1Res(2)- p2Res(2)]*lambda - testRes(2)*testIntensity
%
% Define M = [ [p1Res(1)- p2Res(1)] , -testRes(1) ; [p1Res(2)- p2Res(2)], -testRes(2) ]
% Define b = [-p2Res(1) -p2Res(2)]';
% [lambda testIntensity]' = inv(M)*b;
M = [(p1Res(1)-p2Res(1)),-testRes(1); (p1Res(2)-p2Res(2)),-testRes(2)];
b = [-p2Res(1); -p2Res(2)];
answer = inv(M)*b;
lambda = answer(1);
testIntensity = answer(2);

% Check if the computed scale factors are reasonable
if lambda < 0 || lambda > 1 || testIntensity < 0 || testIntensity > 1
    error('Not possible to compute in range Rayleigh match');
end

% Generate adjusted spds
primaryMixtureSpd = lambda*p1Spd + (1-lambda)*p2Spd;
testAdjustedSpd = testIntensity*testSpd;

% If not monochromatic, add the dark spd
if p.Results.addDarkSpd
    if isempty(p.Results.darkSpd)
        cal = OLGetCalibrationStructure('CalibrationType',...
            getpref('ForcedChoiceCM','currentCal'));
        [~,settingsLength] = size(cal.computed.pr650M);
        darkSpd = OLPrimaryToSpd(cal,zeros(settingsLength,1));
    else
        darkSpd = p.Results.darkSpd;
    end
    primaryMixtureSpd = primaryMixtureSpd+darkSpd;
    testAdjustedSpd = testAdjustedSpd+darkSpd;
end

% Make optional plots of the results
plotResults = false;
if plotResults
    wls = StoWls(p.Results.S);
    primaryLMS = observer.T_cones*primaryMixtureSpd;
    testLMS = observer.T_cones*testAdjustedSpd;
    OLPlotSpdCheck(wls,[primaryMixtureSpd,testAdjustedSpd]);
    legend('Primary','Test');
    title('Adjusted Spds');
    OLPlotConeEffects(primaryLMS,testLMS,'',1);
end
end