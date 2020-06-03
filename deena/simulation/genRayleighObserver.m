function observer = genRayleighObserver(foveal, varargin)
% Generates a simulated observer for a Rayleigh match experiment 
%
% Syntax:
%   genRayleighObserver(); 
%
% Description:
%    Generates a simulated observer based on the Asano model. Designed for
%    use with the OLRayleighMatch program to simulate Rayleigh matching
%    experiments. 
% 
%    If no inputs are given, the program generates an observer with 
%    standard cone fundamentals. However, the user can also enter a vector 
%    of individual difference parameters for the model as a key-value pair
%    (dlens, dmac, dphotopigments, lambdaMaxShifts, and noiseSd). 
%   
%    The output observer struct includes three fields: colorDiffParams 
%    (opponent contrast weighting parameters), coneParams (individual
%    difference/observer parameters), and T_cones (cone fundamentals). 
%
% Inputs:
%     foveal            -Logical indicating whether matches are foveal 
%                       (true) or peripheral (false)   
%
% Outputs:
%     observer          -Struct with observer parameters and cone 
%                        fundamentals. Fields include colorDiffParams, 
%                        coneParams, and T_cones. 
%
% Optional key-value pairs:
%    coneVec           -Numeric vector with cone individual difference
%                       parameters. Includes 9 elements - see 
%                       ObserverVecToParams for a full description. 
 
% History:
%   06/02/20  dce       Wrote initial code
%   06/03/20  dce       Added field size correction


% Parse input 
p = inputParser;
p.addParameter('coneVec', [], @(x) (isvector(x)));
p.parse(varargin{:});

S = [380 2 201];  % Wavelength sampling (following OneLight convention)

observer = struct; 
observer.colorDiffParams = DefaultColorDiffParams('opponentContrast'); 
observer.coneParams = DefaultConeParams('cie_asano'); 
if foveal
    observer.coneParams.fieldSizeDegrees = 2; 
end 
if ~isempty(p.Results.coneVec)
    observer = ObserverVecToParams('basic', p.Results.coneVec, observer);  
end 

observer.T_cones = ComputeObserverFundamentals(observer.coneParams,S); 
end 