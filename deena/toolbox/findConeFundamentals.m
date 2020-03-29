% Finds cone fundamentals for an observer
% Syntax:
%   findConeFundamentals
%
% Description:
%    Returns standard cone fundamentals for an observer. Takes in a 3 x 1
%    vector of lambda max and a 3 x 1 vector of optical density variation,  
%    then produces a matrix of cones that accounts for these individual 
%    differences. Optional key-value pairs allow further refinement of the 
%    model cones.
%
% Inputs: 
%    lambdaMaxes    - 3x1 vector of peak spectral sensitivities
%    dPhotopigments - 3x1 vector of optical densities, expressed as percent
%                     deviation from average/zero
% 
% Outputs:
%    T_cones       - matrix of standard cone fundamentals for the observer
% 
% Optional key-value pairs: 
%    'increment'   - integer distance in nm between wavelength samplings.
%                    Default is 1.
%    'foveal'      - logical indicating whether the observer is making
%                    matches with the fovea (true) or periphery (false). 
%                    Default is false.

% History
%   de  xx/xx/19  - wrote program 
%   de  3/29/20   - added documentation and edits

function T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, varargin)
% Parse input
p = inputParser;
p.addParameter('inc', 1, @(x) (isnumeric(x)));
p.addParameter('foveal', false, @(x) (islogical(x)));
p.parse(varargin{:});

% Observer parameters
if p.Results.foveal
    fieldSizeDegs = 2;
else
    fieldSizeDegs = 10;
end
observerAge = 32;
pupilDiameterMM = 3;

% Individual different parameters for observers
indDiffParams.dlens= 0;
indDiffParams.dmac = 0;
indDiffParams.dphotopigment = dphotopigments;
indDiffParams.lambdaMaxShift = [0 0 0];
indDiffParams.shiftType = 'linear';

S = [380 p.Results.inc ceil(401/p.Results.inc)]; % Wavelength sampling

% Generate the cones
T_cones = EnergyToQuanta(S, ...
    ComputeCIEConeFundamentals(S,fieldSizeDegs,observerAge,...
    pupilDiameterMM,lambdaMaxes, [],[],[],[],[],indDiffParams)')';

% Normalize cones by their maximum 
for ii = 1:size(T_cones,1)
    T_cones(ii,:) = T_cones(ii,:)/max(T_cones(ii,:));
end
end
