function [prob1] = ComputeChoiceLikelihood(params,adaptationLMS,referenceLMS,comparison1LMS,comparison2LMS)
% Compute probability that stim1 will be chosen as closer to reference than stim2
%
% Syntax:
%   [prob1] = ComputeChoiceLikelihood(params,adaptationLMS,referenceLMS,comparison1LMS,comparison2LMS)
%
% Description:
%   This routine simulates a three interval experiment.  The subject's task
%   is to judge whether comparison1 or comparison2 is more similar to the reference.
%
%   All stimuli are specified as spectral power distributions on wavelength
%   sampling given by S.
%
%   The structure params has parameters specifying the observer properties,
%   such as cone spectral sensitivities and color difference model.
%
% Inputs:
%    params                       - Parameter structure describing the observer.
%    adaptationLMS                - Adaptation LMS (3 by 1 vector).
%    referenceLMS                 - Reference LMS (3 by 1 vector)..
%    comparison1LMS               - First LMS (3 by 1 vector)..
%    comparison2LMS               - Second LMS (3 by 1 vector)..
%
% Outputs:
%    prob1                     - The probability that stimulus 1 is judged closer.
%
% Optional key/value pairs
%    None.
%
% See also:
%   ObserverParamsToVec, ObserverVecToParams, ComputeCIEConeFundamentals
%

% History:
%   08/09/19  dhb  Wrote it.
%   08/14/19  dhb  Pass T rather than comnpute locally. It's slow to compute, and this
%                  routine gets called a lot with the same T.

% Compute differences
comparison1Diff = ComputeMatchDiff(params.colorDiffParams,adaptationLMS,referenceLMS,comparison1LMS);
comparison2Diff = ComputeMatchDiff(params.colorDiffParams,adaptationLMS,referenceLMS,comparison2LMS);
diffDiff = comparison2Diff - comparison1Diff;

% Compute likelihood based on differences
prob1 = normcdf(diffDiff,0,params.colorDiffParams.noiseSd);

end