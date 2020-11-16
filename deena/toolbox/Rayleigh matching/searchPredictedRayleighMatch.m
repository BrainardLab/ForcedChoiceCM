function [tSpd,pSpd,tIndex,pIndex,minDiff] =...
    searchPredictedRayleighMatch(testSpds,primarySpds,observer)
% Finds the best possible match for a given set of primary and test lights
%
% Syntax:
%   searchPredictedRayleighMatch(testSpds,primarySpds,observer)
%
% Description:
%    Takes in two matrices of possible primary and test spds used in a 
%    Rayleigh matching simulation, as well as a simulated observer struct. 
%    Conducts an exhaustive search to find which pair of primary and test 
%    spds has the smallest opponent contrast for the observer. Returns the 
%    primary and test spds, as well as their indices in the lights array. 
%
%    Unlike computePredictedRayleighMatch, which computes the ideal match
%    analytically, this routine finds the best match among an available set
%    of lights. This avoids issues of quantization for small adjustment
%    lengths, as well as nonlinearities due to OneLight nominal ->
%    predicted conversion.
%
% Inputs:
%    testSpds     -spdLength x m matrix of possible test spds
%    primarySpds  -spdLength x n matrix of possible primary spds
%    observer     -Structure of a simulated observer (see
%                  genRayleighObserver for details)
%
% Outputs:
%    pSpd        -Predicted primary spd for the match
%    tSpd        -Predicted test spd for the match
%    pIndex      -Index of the selected primary spd in the light settings
%                 array
%    tIndex      -Index of the selected test spd in the light settings
%                 array
%    minDiff     -Opponent contrast norm for the two best lights.
%
% Optional key-value pairs:
%    None

% History
%    dce    08/09/20  - Adapted from findNominalMatch and dhb's
%                       t_AnalyticRayleighMatch
%    dce    08/10/20  - Renamed, changed order of inputs to standardize
%    dce    11/15/20  - Changed to compute opponent contrast of test
%                       relative to primary mixture.

% Initial search settings 
minDiff = Inf;
tIndex = NaN;
pIndex = NaN;

% Calculate cone effects for each primary/test light 
primaryConeEffects = observer.T_cones*primarySpds;
testConeEffects = observer.T_cones*testSpds;

% Loop through lights and find the minimal opponent contrast error (test
% relative to primary mixture)
for primaryTrialInd = 1:size(testSpds,2)
    opponentDiffs = LMSToOpponentContrast(observer.colorDiffParams,...
        primaryConeEffects(:,primaryTrialInd),testConeEffects);
    diffVec = vecnorm(opponentDiffs(1:2,:));
    [diff,testTrialInd] = min(diffVec);
    if (diff < minDiff)
        minDiff = diff;
        tIndex = testTrialInd;
        pIndex = primaryTrialInd;
    end
end 

% Return error if the nominal match is at the limit of available light
% settings. This means a good match may not be possible. 
if (tIndex==size(testSpds,2)) 
    error('Test light is at the limit of available settings');
elseif (pIndex==size(primarySpds,2))
    error('Primary light is at the limit of available settings');
end 

% Find spds of ideal match lights
tSpd = testSpds(:,tIndex);
pSpd = primarySpds(:,pIndex);
testLMS = testConeEffects(:,tIndex);
primaryLMS = primaryConeEffects(:,pIndex);

% Plot cone responses at the minimum error settings
makePlots = false;
if makePlots
    OLPlotConeEffects(primaryLMS,testLMS,'Best Available Match',1);
end
end