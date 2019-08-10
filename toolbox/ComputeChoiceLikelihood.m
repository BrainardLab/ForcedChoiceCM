function [prob1,T] = ComputeChoiceLikelihood(params,S,reference,comparison1,comparison2)
% Compute probability that stim1 will be chosen as closer to reference than stim2
%
% Syntax:
%   [prob1,T] = ComputeChoiceLikelihood(params,S,reference,comparison1,comparison2)
%
% Description:
%   This routine simulates a three interval experiment.  The subject's task
%   is to judge whether comparison1 or comparison2 is more similar to the reference.
%
%   All stimuli are specified as spectral power distributions on wavelength
%   sampling given by S.
%
%   The structure params has parameters specifying the pbserver properties,
%   such as cone spectral sensitivities and color difference model.
%
% Inputs:
%    params                    - Parameter structure describing the observer.
%    S                         - Wavelength sampling for stimuli (PTB convention.)
%    fieldSizeDegrees          - Field size in degrees.
%    ageInYears                - Observer age in years
%    pupilDiameterMM           - Pupil size in mm.
%    reference                 - Reference stimulus spectral power distribution.
%    comparison1               - First comparison spectral power distribution.
%    comparison2               - Second comparison spectral power distribution.
%
% Outputs:
%    prob1                     - The probability that stimulus 1 is judged closer.
%    T                         - Cone fundamentals in energy units
%
% Optional key/value pairs
%    None.
%
% See also:
%   ObserverParamsToVec, ObserverVecToParams, ComputeCIEConeFundamentals
%

% History:
%   08/09/19  dhb  Wrote it.

% Cone fundamentals
T = ComputeObserverFundamentals(params.coneParams,S);

% Compute cone responses
referenceLMS = T*reference;
comparison1LMS = T*comparison1;
comparison2LMS = T*comparison2;

% Compute differences
comparison1Diff = ComputeMatchDiff(params.colorDiffParams,referenceLMS,comparison1LMS);
comparison2Diff = ComputeMatchDiff(params.colorDiffParams,referenceLMS,comparison2LMS);
diffDiff = comparison2Diff - comparison1Diff;

% Compute likelihood based on differences
prob1 = normcdf(diffDiff,0,params.colorDiffParams.noiseSd);

end