function observer = genRayleighObserver(varargin)
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
%    If no key-value pairs are given, the program generates an observer  
%    with standard cone fundamentals, default opponent parameters, and a 
%    2-degree field size. However, the user can also enter a vector of 
%    individual difference parameters as a key-value pair (dlens, dmac, 
%    dphotopigments, and lambdaMaxShifts). The user can also enter a vector
%    of opponent contrast parameters (lumWeight, rgWeight, byWeight, and
%    noise). There is also an option to specify age, field size, and 
%    wavelength sampling. 
%   
%    The output observer struct includes three fields: colorDiffParams 
%    (opponent contrast weighting parameters), coneParams (individual
%    difference/observer parameters), S (wavelength sampling), and 
%    (optionally) T_cones (cone fundamentals). 
% Inputs:
%     none  
%
% Outputs:
%     observer          -Struct with observer parameters and cone 
%                        fundamentals. Fields include colorDiffParams, 
%                        coneParams, S, and (optionally) T_cones. 
%
% Optional key-value pairs:
%    coneVec           -Numeric vector with cone individual difference
%                       parameters. Includes 8 elements - see 
%                       ObserverVecToParams for a full description. 
%    opponentParams    -4-element vector with opponent contrast parameters. 
%                       (1) is the luminance weight, (2) is the rg weight,
%                       (3) is the by weight, and (4) is the noise standard
%                       deviation. Default is [4 2 0.5 0.02].
%    age               -Integer for subject age. Default is 32.
%    fieldSize         -Integer field size in degrees. Default is 2. 
%    calcCones         -Logical indicating whether to calculate cone
%                       fundamentals. Default is true.
%    S                 -Wavelength sampling for cone calculations, in the 
%                       form [start increment numTerms]. Default is 
%                       [380 2 201];  

% History:
%   06/02/20  dce       Wrote initial code
%   06/03/20  dce       Added field size correction
%   06/12/20  dce       Added age and field size as key-value pairs, 
%                       made cone calculation optional 
%   06/29/20  dce       Added S field to observer
%   08/05/20  dce       Edited to have passed-in opponent weights 

% Parse input 
p = inputParser;
p.addParameter('coneVec',zeros(1,8),@(x)(isvector(x)));
p.addParameter('opponentParams',[4 2 0.5 0.02],@(x)(isvector(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x))); 
p.addParameter('calcCones',true,@(x)(islogical(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Input checks 
if length(p.Results.coneVec) ~= 8
    error('Cone parameters vector must have 8 elements'); 
end 
if length(p.Results.opponentParams) ~= 4
    error('Opponent parameters vector must have 4 elements'); 
end 

% Basic setup 
observer = struct; 
observer.colorDiffParams = getColorDiffParams(p.Results.opponentParams); 
observer.coneParams = DefaultConeParams('cie_asano'); 
observer.S = p.Results.S; % Wavelength sampling 

% Set cone individual difference parameters
observer = ObserverVecToParams('basic',[p.Results.coneVec...
    p.Results.opponentParams(4)],observer);

% Reset age if needed
if p.Results.age ~= 32
    observer.coneParams.ageYears = p.Results.age; 
end 

% Reset field size if needed (the default for coneParams is 10 degrees)
if p.Results.fieldSize ~= 10
    observer.coneParams.fieldSizeDegrees =  p.Results.fieldSize; 
end

% Find cone fundamentals
if p.Results.calcCones
    observer.T_cones = ComputeObserverFundamentals(observer.coneParams,...
        p.Results.S);
end 
end 