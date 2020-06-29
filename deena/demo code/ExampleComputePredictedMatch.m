%% Close out stray figures
clear; close all;

%% Set up parameters
%
% Wavelength sampling
S = [400 1 301];
wls = SToWls(S);

% Observer parameters
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
psiParamsStruct.colorDiffParams = DefaultColorDiffParams('opponentContrast');
psiParamsStructRef = psiParamsStruct;
TRef = ComputeObserverFundamentals(psiParamsStructRef.coneParams,S);
T_LM = TRef(1:2,:);

% Test wavelength
testWl = 580;
indexT = find(wls == testWl);
testSpd = zeros(size(wls));
testSpd(indexT) = 1;

% Primaries
primary1Wl = 540;
primary2Wl = 650;
index1 = find(wls == primary1Wl);
index2 = find(wls == primary2Wl);
primary1Spd = zeros(size(wls));
primary1Spd(index1) = 1;
primary2Spd = zeros(size(wls));
primary2Spd(index2) = 1;

% Want this to be true
%   T_LM*(testIntensity)*testSpd == T_LM*(lambda*primary1Spd + (1-lambda)*primary2Spd)
%   0 = T_LM*[(lambda*primary1Spd + (1-lambda)*primary2Spd) - (testIntensity)*testSpd)
%   0 = T_LM*primary1Spd * lambda + T_LM*primary2Spd * (1-lambda) - T_LM*testSpd * testIntensity)
coeff1 = T_LM*primary1Spd;
coeff2 = T_LM*primary2Spd;
coeff3 = T_LM*testSpd;

% So
%    0 = coeff1(1) * lambda + coeff2(1) - coeff2(1)*lambda - coeff3(1)*testIntensity
%    0 = coeff1(2) * lambda + coeff2(2) - coeff2(2)*lambda - coeff3(2)*testIntensity
%
%    -coeff2(1) = [coeff1(1)- coeff2(1]*lambda - coeff3(1)*testIntensity
%    -coeff2(2) = [coeff1(2)- coeff2(2]*lambda - coeff3(2)*testIntensity
%
% Define M = [ [coeff1(1)- coeff2(1)] , -coeff3(1) ; [coeff1(2)- coeff2(2)], -coeff3(2) ]
% Define b = [-coeff2(1) -coeff2(2)]';
% [lambda testIntensity]' = inv(M)*b;
M = [ [coeff1(1)- coeff2(1)] , -coeff3(1) ; [coeff1(2)- coeff2(2)], -coeff3(2) ];
b = [-coeff2(1) -coeff2(2)]';
answer = inv(M)*b;
lambda = answer(1); testIntensity = answer(2);

% Check
primaryMixtureSpd = lambda*primary1Spd + (1-lambda)*primary2Spd;
testAdjustedSpd = testIntensity*testSpd;

primaryLM = T_LM*primaryMixtureSpd;
testLM = T_LM*testAdjustedSpd;
