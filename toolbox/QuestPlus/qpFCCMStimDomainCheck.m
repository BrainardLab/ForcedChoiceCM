function paramsOK = qpFCCMStimDomainCheck(stimParamsVec,type)
% Stimulus parameter check for force-choice color matching
%
% Usage:
%     paramsOK = qpFCCMStimDomainCheck(stimParams,type)
%
% Description:
%     Check whether passed parameters are valid stimuli for forcec choice
%     color matching
%
% Inputs:
%     stimParamsVec     Stimulus parameter vec for the FCCM experiment
%     type              String, parameter type
%
% Output:
%     paramsOK          Boolean, true if parameters are OK and false otherwise.

% 07/22/17  dhb  Wrote it.

%% Assume ok
paramsOK = true;

%% Check that concentration is non-negative
if (psiParams(1) < 0)
    paramsOK = false;
end

%% Check whether boundary parameters are OK
%
% This is signaled by returning NaN when the boundaries are
% not in increasing order.
[boundaries,sortIndex] = sort(psiParams(2:end),'ascend');
nOutcomes = length(boundaries);
if (any(sortIndex ~= 1:nOutcomes))
    paramsOK = false;
end

