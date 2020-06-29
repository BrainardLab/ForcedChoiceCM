function [testAdjustedSpd,primaryMixtureSpd,lambda,testIntensity] =...
    computePredictedRayleighMatch(p1,p2,test,observerParams,varargin)
% Calculates a predicted Rayleigh match analytically

% Syntax:
%   computePredictedRayleighMatch(p1,p2,test,observerParams)
%
% Description:
%    Computes a Rayleigh match analytically when given the associated
%    primary and test wavelengths and the individual difference parameters
%    for the observer (based on the Asano model). Makes a monochromatic spd
%    for each of the wavelengths, then solves a matrix equation to find the
%    primary ratio and test intensity which would lead to a match. There is
%    also an option to plot the chosen match spds and their cone effects.
%
% Inputs:
%    p1              -Integer wavelength of the first primary light (nm).
%    p2              -Integer wavelength of the second primary light (nm).
%    test            -Integer wavelength of the test light (nm).
%    observerParams  -Nine-element vector with the observer's individual
%                      difference parameters (see findObserverParameters 
%                      for a full description)
%                
% Outputs:
%    primaryMixtureSpd    -Predicted primary spd for the match
%    testAdjustedSpd      -Predicted test spd for the match
%    lambda               -Scale factor for primary 1 proportion, between 0
%                          and 1.
%    testIntensity        -Scale factor for test intensity, between 0 and 1
%
% Optional key-value pairs:
%    age               -Integer for subject age. Default is 32.
%    fieldSize         -Integer field size in degrees. Default is 2. 

% History
%    dce    6/29/20   -Adapted from example code from dhb 

%% Parse input 
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.parse(varargin{:});

%% Set up parameters
% Wavelength sampling
S = [400 1 301];
wls = SToWls(S);

% Generate monochromatic spds 
p1Spd = zeros(size(wls)); 
p1Spd(wls==p1) = 1; 

p2Spd = zeros(size(wls)); 
p2Spd(wls==p2) = 1; 

testSpd = zeros(size(wls));
testSpd(wls==test) = 1; 

% Construct a simulated observer with the specified parameters 
observer = genRayleighObserver('fieldSize',p.Results.fieldSize,'age',...
    p.Results.age,'coneVec', observerParams,'S',S);
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
    error('Not possible to compute Rayleigh match for the given primaries'); 
end 

% Generate adjusted spds 
primaryMixtureSpd = lambda*p1Spd + (1-lambda)*p2Spd;
testAdjustedSpd = testIntensity*testSpd;
primaryLMS = observer.T_cones*primaryMixtureSpd;
testLMS = observer.T_cones*testAdjustedSpd;

% Make optional plots of the results 
plotResults = false;
if plotResults 
    OLPlotSpdCheck(wls,[primaryMixtureSpd,testAdjustedSpd]);
    legend('Primary','Test'); 
    title('Adjusted Spds'); 
    OLPlotConeEffects(primaryLMS,testLMS,'',1); 
end 
end 