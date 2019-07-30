function prob1 = ComputeChoiceLikelihood(params,S,reference,comparison1,comparison2)
% Compute probability that stim1 will be chosen as closer to reference than stim2
%
% Syntax:
%   prob1 = ComputeChoiceLikelihood(params,S,reference,comparison1,comparison2)
%
% Description:
%   This routine simulates a three interval experiment.  The subject's task
%   is to judge whether comparison1 or comparison2 is more similar to the reference.
%
%   All stimuli are specified as spectral power distributions on wavelength
%   sampling given by S.
%
%   The structure params specifies the cone spectral sensitivities of the
%   observer, as well as a noise model.
%
% Inputs:
%    params                    - Parameter structure describing the
%                                observer.
%    S                         - Wavelength sampling for stimuli (PTB
%                                convention.)
%    reference                 - Reference stimulus spectral power
%                                distribution.
%    comparison1               - First comparison spectral power
%                                distribution.
%    comparison2               - Second comparison spectral power
%                                distribution.
%
% Outputs:
%    prob1                     - The probability that stimulus 1 is judged
%                                closer.
%
% Optional key/value pairs
%    None.

% Check inputs and parse

% Get cone spectral sensitivities

% Compute cone responses

% Convert to cone contrast space

% Compute likelihood based on distance

end