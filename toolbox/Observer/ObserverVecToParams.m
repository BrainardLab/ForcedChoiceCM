function params = ObserverVecToParams(x)
% Convert vector of observer parameters to a structure
%
% Synopsis:
%   params = ObserverVecToParams([x])
%
% Description:
%   Our goal is to used forced choice color similarity judgments to
%   determine observer parameters.  Sometimes we want those parameters
%   as a vector, and sometimes in a structure.  This routine goes from
%   the vector to the structure.
%
%   This illustrates the transformation.
%       params.indDiffParams.dlens = x(1);
%       params.indDiffParams.dmac = x(2);
%       params.indDiffParams.dphotopigment(1) = x(3);
%       params.indDiffParams.dphotopigment(2) = x(4);
%       params.indDiffParams.dphotopigment(3) = x(5);
%       params.indDiffParams.lambdaMaxShift(1) = x(6);
%       params.indDiffParams.lambdaMaxShift(2) = x(7);
%       params.indDiffParams.lambdaMaxShift(3) = x(8);
%
%
% Inputs:
%   x                       - Parameters as vector.  If not passed, it is
%                             set to all zeros.
%
% Outputs:
%   params                  - Parameter structure.
%
% Optional key value pairs:
%   None.
%
% See also:
%   ObserverParamsToVec, ComputeCIEConeFundamentals
%

% History:
%   08/09/19  dhb  Wrote it, because I have to do one fun thing this summer.

% Examples:
%{
    x = (1:8)';
    params = ObserverVecToParams(x);
    params.indDiffParams
    x1 = ObserverParamsToVec(params)
    if (any(x - x1) ~= 0)
        error('Routines do not properly self invert');
    end
%}

% This routine might get called a lot, so don't use the
% wonderful but lugubrious input parser.
if (nargin == 0 | isempty(x))
    x = zeros(8,1);
end

params.indDiffParams.dlens = x(1);
params.indDiffParams.dmac = x(2);
params.indDiffParams.dphotopigment(1) = x(3);
params.indDiffParams.dphotopigment(2) = x(4);
params.indDiffParams.dphotopigment(3) = x(5);
params.indDiffParams.lambdaMaxShift(1) = x(6);
params.indDiffParams.lambdaMaxShift(2) = x(7);
params.indDiffParams.lambdaMaxShift(3) = x(8);
params.indDiffParams.shiftType = 'linear';
