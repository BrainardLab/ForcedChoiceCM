% Returns cone fundamentals for an observer. Takes in a 3 x 1 vector of
% lambda max and a 3 x 1 vector of optical density variation. Optional
% key-value pairs include'increment', an integer for the distance in nm
% between wavelength samplings (default is 1); and 'foveal', which 
% indicates whether the observer is making matches with the fovea (true) or 
% the periphery (false) - (default is false).  

function T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, varargin)
% Parse input
p = inputParser;
p.addParameter('inc', 1, @(x) (isnumeric(x)));
p.addParameter('foveal', false, @(x) (islogical(x))); 
p.parse(varargin{:});
inc = p.Results.inc; 

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

% % Generate the cones 
S = [380 inc ceil(401/inc)]; % Wavelength sampling 
T_cones = EnergyToQuanta(S, ...
    ComputeCIEConeFundamentals(S,fieldSizeDegs,observerAge,...
    pupilDiameterMM,lambdaMaxes, [],[],[],[],[],indDiffParams)')';
for ii = 1:size(T_cones,1)
    T_cones(ii,:) = T_cones(ii,:)/max(T_cones(ii,:));
end
end 
