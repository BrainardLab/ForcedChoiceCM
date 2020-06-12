function observer = genRayleighObserver(varargin)
% Generates a simulated observer for a Rayleigh match experiment 
%
% Syntax:
%   genRayleighObserver(foveal); 
%
% Description:
%    Generates a simulated observer based on the Asano model. Designed for
%    use with the OLRayleighMatch program to simulate Rayleigh matching
%    experiments. 
% 
%    If no key-alue pairs are given, the program generates an observer with 
%    standard cone fundamentals and a field size specified by the "foveal" 
%    input. However, the user can also enter a vector of individual 
%    difference parameters as a key-value pair (dlens, dmac, 
%    dphotopigments, lambdaMaxShifts, and noiseSd). 
%   
%    The output observer struct includes three fields: colorDiffParams 
%    (opponent contrast weighting parameters), coneParams (individual
%    difference/observer parameters), and T_cones (cone fundamentals). 
%
% Inputs:
%     none  
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
%    age               -Integer for subject age. Default is 32.
%    fieldSize         -Integer field size in degrees. Default is 10. 
%    calcCones         -Logical indicating whether to calculate cone
%                       fundamentals. Default is true.

% History:
%   06/02/20  dce       Wrote initial code
%   06/03/20  dce       Added field size correction
%   06/12/20  dce       Added age and field size as key-value pairs, 
%                       made cone calculation optional 


% Parse input 
p = inputParser;
p.addParameter('coneVec',[],@(x)(isvector(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',10,@(x)(isnumeric(x))); 
p.addParameter('calcCones',true,@(x)(islogical(x)));
p.parse(varargin{:});

S = [380 2 201];  % Wavelength sampling (following OneLight convention)

observer = struct; 
observer.colorDiffParams = DefaultColorDiffParams('opponentContrast'); 
observer.coneParams = DefaultConeParams('cie_asano'); 

% Reset age if needed
if p.Results.age ~= 32
    observer.coneParams.ageYears = p.Results.age; 
end 

% Reset field size if needed (the default is 10 degrees)
if p.Results.fieldSize ~= 10
    observer.coneParams.fieldSizeDegrees =  p.Results.fieldSize; 
end

% Set individual difference parameters
if ~isempty(p.Results.coneVec)
    observer = ObserverVecToParams('basic', p.Results.coneVec, observer);  
end 

% Cone fundamentals
if p.Results.calcCones
    observer.T_cones = ComputeObserverFundamentals(observer.coneParams,S);
end 
end 