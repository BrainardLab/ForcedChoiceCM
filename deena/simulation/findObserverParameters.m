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
%    constraints can be entered as key-value pairs (see below). Also has an
%    option to plot the optimal cone response to the test and match lights.
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
%    'restrictBySD'  -Logical. If true, adds lower and upper bounds on all
%                     paramters to keep them within three standard
%                     deviations of their means. Default is true. 
%    'S'             -Wavelength sampling for cone calculations, in the
%                     form [start delta nTerms]. Default is [380 2 201];
%    'errScalar'     -Integer for scaling the match error, to improve search.
%                     Default is 100.
%    'initialConeParams' -1x8 numerical vector of cone individual  
%                         difference parameters. Default is zeros(1,8);
%    'opponentParams'    -1x4 numerical vector of opponent contrast
%                         parameters. Default is [4 2 0.5 0.02].

% History:
%   06/12/20  dce       Wrote it.
%   07/06/20  dce       Added S setting and appropriate length checks.
%   07/21/20  dce       Added option to set dlens to 0.
%   07/24/20  dce       Added option to set dmac to 0.
%   08/05/20  dce       Added opponent contrast params 

%% Initial Setup
% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',false,@(x)(islogical(x)));
p.addParameter('dmac0',false,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('initialConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[4 2 0.5 0.02],@(x)(isvector(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Input checks
[spdLength,~] = size(testSpds);
if length(SToWls(p.Results.S)) ~= spdLength
    error('Chosen S value is incompatible with spd length');
end

if p.Results.LMEqualOD && p.Results.dlens0
    error('Only one equality constraint can be selected');
end 

% Generate a standard observer with the given initial values
observer = genRayleighObserver('fieldSize',p.Results.fieldSize,'age',...
    p.Results.age,'calcCones',false,'coneVec', p.Results.initialConeParams,...
    'S',p.Results.S,'opponentParams',p.Results.opponentParams);

%% Restrictions on parameters
Aeq = [];              % Equality constraint
Beq = [];              % Equality constraint 
lb = -Inf*ones(1,8);   % Lower bounds
ub = Inf*ones(1,8);    % Upper bounds

% Set lower and upper bounds based on the standard deviations for the 8
% parameters (from Asano 2015). The standard deviations are expressed as
% percent deviations from the mean, except for last three parameters
% (lambda max shifts) which are expressed as deviations in nm.
sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3]; % Standard deviations
scaleFactor = 3;    % Set limits at 3 standard deviations from the mean
if p.Results.restrictBySd
    lb = -1*scaleFactor*sds;
    ub = scaleFactor*sds;
end

% Constrain lens density variation to 0 
if p.Results.dlens0
    lb(1) = p.Results.initialConeParams(1);
    ub(1) = p.Results.initialConeParams(1);
end

% Constrain lens density variation to 0 
if p.Results.dmac0
    lb(2) = p.Results.initialConeParams(2);
    ub(2) = p.Results.initialConeParams(2);
end

% Constrain L and M OD to be equal
if p.Results.LMEqualOD
    Aeq = [0 0 1 -1 0 0 0 0];
    Beq = 0;
end

%% Optimization procedure
% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter',...
    'LargeScale','off','Algorithm','active-set');

% Find optimal parameters
params = fmincon(@(x)findMatchError(x,observer,testSpds,primarySpds,...
    'S',p.Results.S,'errScalar',p.Results.errScalar),....
    p.Results.initialConeParams(1:8),[],[],Aeq,Beq,lb,ub,[],options);

% Reset observer with the final parameters
observer = ObserverVecToParams('basic', ...
    [params observer.colorDiffParams.noiseSd],observer);
observer.T_cones = ComputeObserverFundamentals(observer.coneParams,...
    p.Results.S);

% What is the error?
error = findMatchError(params,observer,testSpds,primarySpds,'S',...
    p.Results.S,'errScalar',p.Results.errScalar)/p.Results.errScalar;

%% Optional - plot LMS response
plotLMS = false;
if plotLMS
    for i = 1:nMatches
        testCones = observer.T_cones*testSpds(:,i);
        primaryCones = observer.T_cones*primarySpds(:,i);
        cones = [testCones(1), primaryCones(1); testCones(2),...
            primaryCones(2); testCones(3), primaryCones(3)];
        figure();
        bar(cones);
        hold on;
        names ={'L'; 'M'; 'S'};
        set(gca,'xticklabel', names)
        ylabel('Relative Response Intensity');
        legend('Test Response', 'Primary Response');
        theTitle = sprintf('Cone Responses for Optimal Parameters, Match %g', i);
        title(theTitle);
    end
end
end