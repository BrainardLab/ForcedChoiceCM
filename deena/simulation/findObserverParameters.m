function [params, error, observer] = findObserverParameters(testSpd,primarySpd,varargin)
% Finds optimal cone paramters for an observer based on a Rayleigh match
% they made

% Syntax:
%   findObserverParams(testSpd,primarySpd)
%
% Description:
%    Takes in a pair of test and primary spds which a user identified as a
%    Rayleigh match. Then, finds the Asano individual difference parameters
%    which lead these two lights to cause the most similar cone responses.
%    Does this by using fmincon to minimize the opponent contrast of the
%    test and match lights. Various optimization constraints can be entered
%    as key-value pairs (see below). Also creates optional plots of the
%    optimal cone response to the test and match lights.
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
%    testSpd     -201x1 vector representation of the predicted spd for the
%                 chosen test light
%    primarySpd  -201x1 vector representation of the predicted spd for the
%                 chosen primary light
%
% Outputs:
%    params      -1x8 vector containing optimized individual difference
%                 paramters.
%    error       -Vector length of the luminance and RG opponent contrast
%                 terms, when optimal cone paramters are used.
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
%    'restrictBySD'  -Logical. If true, adds lower and upper bounds on all
%                     paramters to keep them within three standard
%                     deviations of their means. Default is false.

% History:
%   06/12/20  dce       Wrote it.

% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('restrictBySD',false,@(x)(islogical(x)));
p.parse(varargin{:});

% Generate a standard observer, and set all initial individual difference
% parameters to 0. (this is the default for the standard observer).
observer = genRayleighObserver('fieldSize',p.Results.fieldSize,...
    'age',p.Results.age,'calcCones',false);
obsParams = ObserverParamsToVec('basic',observer);
initialParams = obsParams(1:8);

% Lower and upper bounds, other restrictions
Aeq = [];
Beq = [];
lb = [];
ub = [];
if p.Results.LMEqualOD
    Aeq = [0 0 1 -1 0 0 0 0];
    Beq = 0;
end

% Set fmincon options
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','off',...
    'LargeScale','off','Algorithm','active-set');

% Find optimal parameters
params = fmincon(@(x)findMatchError(x,observer,testSpd,primarySpd),...
    initialParams,[],[],Aeq,Beq,lb,ub,[],options);

% Reset observer with the final parameters
observer = ObserverVecToParams('basic', ...
    [params observer.colorDiffParams.noiseSd], observer);
observer.T_cones = ComputeObserverFundamentals(observer.coneParams,...
    [380 2 201]);

% What is the error?
error = findMatchError(params,observer,testSpd,primarySpd);

% Plot LMS response
plotLMS = true;
if plotLMS
    testCones = observer.T_cones*testSpd;
    primaryCones = observer.T_cones*primarySpd;
    cones = [testCones(1), primaryCones(1); testCones(2),...
        primaryCones(2); testCones(3), primaryCones(3)];
    figure();
    bar(cones);
    hold on;
    names ={'L'; 'M'; 'S' };
    set(gca,'xticklabel', names)
    ylabel('Relative Response Intensity');
    legend('Test Response', 'Primary Response');
    title('Cone Responses for Optimal Parameters');
end
end