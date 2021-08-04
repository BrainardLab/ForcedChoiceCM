function [params,error,observer] = findObserverParameters(testSpds,primarySpds,varargin)
% Finds optimal cone paramters for an observer based on a Rayleigh match
% they made

% Syntax:
%   findObserverParameters(testSpds,primarySpds)
%
% Description:
%    Takes in an array of test and primary spd pairs which a user
%    identified as a Rayleigh match. Then, finds the Asano individual
%    difference parameters which cause these light pairs to have the most
%    similar cone responses. Does this by using fmincon to minimize the
%    opponent contrast of the test and match lights. Various optimization
%    constraints can be entered as key-value pairs (see below), or limits
%    can be entered directly.
%
%    The individual difference paramters are as follows:
%       observer.coneParams.indDiffParams.dlens = x(1);
%       observer.coneParams.indDiffParams.dmac = x(2);
%       observer.coneParams.indDiffParams.dphotopigment(1) = x(3);
%       observer.coneParams.indDiffParams.dphotopigment(2) = x(4);
%       observer.coneParams.indDiffParams.dphotopigment(3) = x(5);
%       observer.coneParams.indDiffParams.lambdaMaxShift(1) = x(6);
%       observer.coneParams.indDiffParams.lambdaMaxShift(2) = x(7);
%       observer.coneParams.indDiffParams.lambdaMaxShift(3) = x(8);
%
% Inputs:
%    testSpds     -Vector representation of the predicted spds for
%                  the chosen test lights
%    primarySpds  -Vector representation of the predicted spds for
%                  the chosen primary lights
%
% Outputs:
%    params      -1x8 vector containing optimized individual difference
%                 paramters.
%    error       -Vector length of the luminance and RG opponent contrast
%                 terms when the final cone paramters are used.
%    observer    -Simulated observer struct which includes the optimized
%                 individual difference parameters. Follows the format of
%                 genRayleighObserver().
%
% Optional key-value pairs:
%    'age'           -Observer age in years. Default is 32.
%    'fieldSize'     -Field size in degrees. Default is 2.
%    'LMEqualOD'     -Logical. If true, constrains the L and M cone optical
%                     densities (params(3) and params(4)) to be equal.
%                     Default is false.
%    'dlens0'        -Logical. If true, constrains the lens pigment density
%                     (params(1)) to be 0. Default is false.
%    'dmac0'         -Logical. If true, constrains the macular pigment 
%                     density(params(2)) to be 0. Default is false.
%    'OD0'           -Logical. If true, constrains the optical
%                     density(params(3:5)) to be 0. Default is false.
%    'lambdaMax0'    -Logical. If true, constrains lambda max
%                     (params(6:8)) to be 0. Default is false.
%    'S0'            -Logical. If true, constrains S lambda max and optical
%                     density (params 5 and 8) to be 0. Default is true.
%    'restrictBySD'  -Logical. If true, adds lower and upper bounds on all
%                     paramters to keep them within a specified number of 
%                     standard deviations of their means. Default is true. 
%    'S'             -Wavelength sampling for cone calculations, in the
%                     form [start delta nTerms]. Default is [380 2 201];
%    'errScalar'     -Integer for scaling the match error, to improve search.
%                     Default is 100.
%    'initialConeParams' -1x8 numerical vector of cone individual  
%                         difference parameters. Default is zeros(1,8);
%    'opponentParams'    -1x4 numerical vector of opponent contrast
%                         parameters. Default is [40.3908 205.7353 62.9590 1.0000].
%    'minimizeConeErr'   -Logical. If true, minimizes cone exictation error
%                         instead of opponent contrast difference. Default 
%                         is false.
%    'lowerBounds'   -1x8 vector with lower bounds. Overwrites key-value
%                     preferences. Default is [].
%    'upperBounds'   -1x8 vector with upper bounds. Overwrites key-value
%                     preferences. Default is [].
%    'AEq'           -1x8 vector with linear equality constraints. 
%                     Overwrites key-value preferences. Default is [].
%    'BEq'           -1x8 vector with linear equality sum. Overwrites 
%                     key-value preferences. Default is [].
%    'sdDensity'     -Number of allowed standard deviations for density 
%                     parameters (1:5) in fit. Default is 3.
%    'sdLambdaMax'   -Number of allowed standard deviations for lambda max
%                     parameters (6:8) in fit. Default is 3.
%    'matchErrFun'   -Handle to function that is used to evaluate the fit
%                     error in the minimization.

% History:
%   06/12/20  dce       Wrote it.
%   07/06/20  dce       Added S setting and appropriate length checks.
%   07/21/20  dce       Added option to set dlens to 0.
%   07/24/20  dce       Added option to set dmac to 0.
%   08/05/20  dce       Added opponent contrast params 
%   11/01/20  dce       Restricted to within 2 standard deviations
%   01/01/21  dce       Restricted to within 3 standard deviations
%   03/01/21  dce       Added option to lock L and M optical densities at
%                       0.
%   04/30/21  dce       Added S locking as a key value pair.
%   05/09/21  dce       Added option to minimize cone excitation difference
%                       instead of opponent contrast.
%   06/08/21  dce       Renamed parameters, added option to constrain
%                       lambda max to 0
%   06/30/21  dce       Added option to directly input lower and upper
%                       bounds and equality constraints - useful for cross 
%                       validation. Note that these settings overwrite
%                       key-value pair results.
%   07/05/21  dce       Added number of sds as key-value pairs

%% Initial Setup
% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',false,@(x)(islogical(x)));
p.addParameter('dmac0',false,@(x)(islogical(x)));
p.addParameter('OD0',false,@(x)(islogical(x)));
p.addParameter('lambdaMax0',false,@(x)(islogical(x)));
p.addParameter('S0',true,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('initialConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.addParameter('minimizeConeErr',false,@(x)(islogical(x)));
p.addParameter('lowerBounds',[],@(x)(isnumeric(x)));
p.addParameter('upperBounds',[],@(x)(isnumeric(x)));
p.addParameter('AEq',[],@(x)(isnumeric(x)));
p.addParameter('BEq',[],@(x)(isnumeric(x)));
p.addParameter('sdDensity',3,@(x)(isnumeric(x)));
p.addParameter('sdLambdaMax',3,@(x)(isnumeric(x)));
p.addParameter('matchErrorFun',@findMatchError,@(x) isa(x,'function_handle'));
p.parse(varargin{:});

% Input checks
spdLength = size(testSpds,1);
if length(SToWls(p.Results.S)) ~= spdLength
    error('Chosen S value is incompatible with spd length');
end
% Check if lower bounds were provided without upper bounds, or if partial
% equality constraints were inputted
if (isempty(p.Results.lowerBounds) || isempty(p.Results.upperBounds)) && ...
        ~(isempty(p.Results.lowerBounds) && isempty(p.Results.upperBounds))
    error('Both lower and upper bounds must be provided');
end 
if (isempty(p.Results.AEq) || isempty(p.Results.BEq)) && ...
        ~(isempty(p.Results.AEq) && isempty(p.Results.BEq))
    error('Both AEq and BEq must be provided');
end 

% Generate a standard observer with the given initial values
observer = genRayleighObserver('fieldSize',p.Results.fieldSize,'age',...
    p.Results.age,'calcCones',false,'coneVec', p.Results.initialConeParams,...
    'S',p.Results.S,'opponentParams',[40.3908 205.7353 62.9590 1.0000]);

%% Restrictions on parameters
Aeq = [];              % Equality constraint
Beq = [];              % Equality constraint 
lb = -Inf*ones(1,8);   % Lower bounds
ub = Inf*ones(1,8);    % Upper bounds

% Set lower and upper bounds, either with explicit input or using key-value
% pairs for specific constraints
if ~isempty(p.Results.lowerBounds)  % Explicit input of bounds
    lb = p.Results.lowerBounds;
    ub = p.Results.upperBounds;
else       % Key-value pairs with specified constraints
    % Set lower and upper bounds based on the standard deviations for the 8
    % parameters (from Asano 2015). The standard deviations are expressed as
    % percent deviations from the mean, except for last three parameters
    % (lambda max shifts) which are expressed as deviations in nm.
    sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3]; % Standard deviations
    
    % Set limits at the specified number of standard deviations from the mean
    scaleFactors = [repmat(p.Results.sdDensity,1,5)...
        repmat(p.Results.sdLambdaMax,1,3)];  
    if p.Results.restrictBySd
        lb = -1.*scaleFactors.*sds;
        ub = scaleFactors.*sds;
    end
    
    % Constrain lens density variation to 0
    if p.Results.dlens0
        lb(1) = p.Results.initialConeParams(1);
        ub(1) = p.Results.initialConeParams(1);
    end
    
    % Constrain macular pigment density variation to 0
    if p.Results.dmac0
        lb(2) = p.Results.initialConeParams(2);
        ub(2) = p.Results.initialConeParams(2);
    end
    
    % Constrain optical density variation to 0
    if p.Results.OD0
        lb(3) = p.Results.initialConeParams(3);
        ub(3) = p.Results.initialConeParams(3);
        lb(4) = p.Results.initialConeParams(4);
        ub(4) = p.Results.initialConeParams(4);
        lb(5) = p.Results.initialConeParams(5);
        ub(5) = p.Results.initialConeParams(5);
    end
    
    % Constrain lambda max variation to 0
    if p.Results.lambdaMax0
        lb(6) = p.Results.initialConeParams(6);
        ub(6) = p.Results.initialConeParams(6);
        lb(7) = p.Results.initialConeParams(7);
        ub(7) = p.Results.initialConeParams(7);
        lb(8) = p.Results.initialConeParams(8);
        ub(8) = p.Results.initialConeParams(8);
    end
    
    % Constrain S cone parameters to 0
    if p.Results.S0
        lb(5) = p.Results.initialConeParams(5);
        ub(5) = p.Results.initialConeParams(5);
        lb(8) = p.Results.initialConeParams(8);
        ub(8) = p.Results.initialConeParams(8);
    end
end

% Set equality constraints, either through explicit input or through
% key-value pairs
if ~isempty(p.Results.lowerBounds) % Explicit input has been provided
    Aeq = p.Results.AEq;
    Beq = p.Results.BEq;
elseif p.Results.LMEqualOD
    Aeq = [0 0 1 -1 0 0 0 0];
    Beq = 0; 
end

%% Optimization procedure
% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter',...
    'LargeScale','off','Algorithm','active-set');

% Find optimal parameters
params = fmincon(@(x)p.Results.matchErrorFun(x,observer,testSpds,primarySpds,...
    'S',p.Results.S,'errScalar',p.Results.errScalar,'findConeErr',...
    p.Results.minimizeConeErr),p.Results.initialConeParams(1:8),...
    [],[],Aeq,Beq,lb,ub,[],options);

% Reset observer with the final parameters
observer = ObserverVecToParams('basic', ...
    [params observer.colorDiffParams.noiseSd],observer);
observer.T_cones = ComputeObserverFundamentals(observer.coneParams,...
    p.Results.S);

% What is the error?
error = p.Results.matchErrorFun(params,observer,testSpds,primarySpds,'S',...
    p.Results.S,'errScalar',p.Results.errScalar,'findConeErr',...
    p.Results.minimizeConeErr)/p.Results.errScalar;
end