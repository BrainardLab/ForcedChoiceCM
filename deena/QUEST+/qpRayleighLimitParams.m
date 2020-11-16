function paramsOK = qpRayleighLimitParams(stimParams,paramLimits)
% qpRayleighLimitParams  Parameter check for qpRayleighSim
%
% Usage:
%     paramsOK = qpRayleighLimitParams(psiParams,paramLimits)
%
% Description:
%     Removes combinations of lambda, test intensity, and test wavelength
%     which fall outside the matching range for normal observers. This
%     allows the QUEST+ Rayleigh mtaching routine to be more efficient.
%
%     Takes in paramLimits, a nTestWls x 5 matrix. Each row represents a
%     given test wavelength and the limits which are associated with it. 
%     The columns are arranged as follows: 
%     [test wl, min lambda, max lambda, min test intensity, max test
%     intensity]
%
% Inputs:
%     stimParams      Matrix of possible stimulus combinations. See 
%                     qpRayleighSim.
%     paramLimits     See description above.
%
% Output:
%     paramsOK       Boolean, true if parameters are OK and false otherwise.

% 11/13/20  dce  Wrote it, adapted from qpPFCircularParamsCheck.

%% Assume ok
paramsOK = true;

%% Find which wavelength we are using
testWl = stimParams(3);
limitsInd = find(paramLimits(:,1)==testWl,1);
if isempty(limitsInd)
    error('Limits not provided for specifed psiParams wavelength %g',testWl);
end 

%% Check whether the parameter falls within limits
if stimParams(1) < paramLimits(limitsInd,2) || ...
        stimParams(1) > paramLimits(limitsInd,3) || ...
        stimParams(2) < paramLimits(limitsInd,4) || ...
        stimParams(2) > paramLimits(limitsInd,5)
    paramsOK = false;
end 
end 