% Script which demonstrates how to run a simple Rayleigh matching
% experiment using QUEST+. Demonstrates the use of qpRayleighSim with
% default parameters. 

% History
%   10/27/20  dce       Wrote it 

% Define which wavelengths we will use in the Rayleigh matching experiment. 
%   p1Wl and p2Wl are the peaks of the two primaries in the primary mixture. 
%   testWls is a set of peak wavelengths for possible test lights.
p1Wl = 670;
p2Wl = 560;
testWls = [570 590 610 630 650];

% Simulated cone parameters for the observer. This is a vector of eight 
% individual difference parameters: lens pigment density, macular pigment 
% density, L/M/S photopigment densities, and L/M/S peak spectral
% sensitivities (lambda max) 
simConeParams = [0 0 20 10 0 -3 2 0];
baseConeParams = zeros(1,8);           % Baseline parameters
coneParamsToVary = [0 0 1 1 0 1 1 0];  % The parameters we allow QUEST+ to vary


% Additional parameters
noiseScaleFactor = 3;   % Scale factor for observer noise
lambdaRef = 0.8;        % Primary ratio used as a reference for opponent contrsat calculations
S = [400 1 301];        % Wavelength sampling
nTrials = 20;           % Number of QUEST+ trials
nObservers = 1;         % Number of observers

% Run the simulation
[questData,psiParams,psiFit] = qpRayleighSim('test',nObservers,nTrials,...
    baseConeParams,coneParamsToVary,noiseScaleFactor,p1Wl,p2Wl,testWls,...
    'precomputeQuest',true,'lambdaRef',lambdaRef,'sampledObservers',simConeParams);
