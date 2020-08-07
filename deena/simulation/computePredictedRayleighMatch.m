function [testAdjustedSpd,primaryMixtureSpd,testIntensity,lambda] =...
    computePredictedRayleighMatch(p1Spd,p2Spd,testSpd,observerParams,...
    opponentParams,varargin)
% Calculates a predicted Rayleigh match analytically
%
% Syntax:
%   computePredictedRayleighMatch(p1,p2,test,observerParams,opponentParams)
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
%
% Inputs:
%    p1Spd           -nx1 spd of the first primary light.
%    p2Spd           -nx1 of the second primary light.
%    testSpd         -nx1 of the test light.
%    observerParams  -Eight-element vector with the observer's individual
%                     difference parameters (see findObserverParameters
%                     for a full description)
%    opponentParams  -4-element vector with opponent contrast parameters. 
%                     (1) is the luminance weight, (2) is the rg weight,
%                     (3) is the by weight, and (4) is the noise standard
%                     deviation. 
%
% Outputs:
%    primaryMixtureSpd    -Predicted primary spd for the match
%    testAdjustedSpd      -Predicted test spd for the match
%    lambda               -Scale factor for primary 1 proportion, between 0
%                          and 1.
%    testIntensity        -Scale factor for test intensity, between 0 and 1
%
% Optional key-value pairs:
%    'age'             -Integer for subject age. Default is 32.
%    'fieldSize'       -Integer field size in degrees. Default is 2.
%    'monochromatic'   -Logical indicating that the passed spds are
%    'darkSpd'         -nx1 dark spd vector. If not specified, default is
%                       [].
%    'S'               -Wavelength sampling for cone calculations, in the
%                       form [start delta nTerms]. Default is [380 2 201],
%                       which is the OneLight convention. Note that S can
%                       only vary from the OneLight convention when
%                       monochromatic lights are used. Otherwise, the
%                       program will throw an error.

% History
%    dce    6/29/20   -Adapted from example code from dhb
%    dce    7/01/20   -Added noise
%    dce    7/03/20   -Add noise to lambda and testIntensity, not cone res
%    dce    7/06/20   -Added option to use OL spds, added scale factors,
%                      made no noise the default
%    dce    7/14/20   -Added predicted matching option
%    dce    7/16/20   -Separated spd creation into a separate file
%    dce    7/22/20   -Added option to pass in dark spd

%% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('monochromatic',false,@(x)(islogical(x)));
p.addParameter('darkSpd',[],@(x)(isnumeric(x)));
p.parse(varargin{:});

%% Set up parameters
% Construct a simulated observer with the specified parameters
observer = genRayleighObserver('fieldSize',p.Results.fieldSize,'age',...
    p.Results.age,'coneVec',observerParams,'opponentParams',opponentParams,...
    'S',p.Results.S);
T_LM = observer.T_cones(1:2,:);         % L and M cone fundamentals

% Cone responses to the three base spds. Each vector contains two elements
% - one for the L cone, one for the M cone.
p1Res = T_LM*p1Spd;
p2Res = T_LM*p2Spd;
testRes = T_LM*testSpd;

%% Set up a matrix equation to solve for the optimal primary ratio and
%% test intensity
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
    error('Not possible to compute Rayleigh match');
end

% Generate adjusted spds
primaryMixtureSpd = lambda*p1Spd + (1-lambda)*p2Spd;
testAdjustedSpd = testIntensity*testSpd;

% If not monochromatic, add the dark spd
if ~p.Results.monochromatic
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
    primaryLMS = observer.T_cones*primaryMixtureSpd;
    testLMS = observer.T_cones*testAdjustedSpd;
    OLPlotSpdCheck(wls,[primaryMixtureSpd,testAdjustedSpd]);
    legend('Primary','Test');
    title('Adjusted Spds');
    OLPlotConeEffects(primaryLMS,testLMS,'',1);
end
end