function x = ObserverParamsToVec(params)
% Convert structure of observer parameters to vector form
%
% Synopsis:
%   params = ObserverVecToParams(x)
%
% Description:
%   Our goal is to used forced choice color similarity judgments to
%   determine observer parameters.  Sometimes we want those parameters
%   as a vector, and sometimes in a structure.  This routine goes from
%   the structure to the vector.
%
%   This illustrates the transformation, in the reverse direction
%       params.indDiffParams.dlens = x(1);
%       params.indDiffParams.dmac = x(2);
%       params.indDiffParams.dphotopigment(1) = x(3);
%       params.indDiffParams.dphotopigment(2) = x(4);
%       params.indDiffParams.dphotopigment(3) = x(5);
%       params.indDiffParams.lambdaMaxShift(1) = x(6);
%       params.indDiffParams.lambdaMaxShift(2) = x(7);
%       params.indDiffParams.lambdaMaxShift(3) = x(8);
%
% Inputs:
%   params                  - Parameter structure.
%
% Outputs:
%   x                       - Parameters as vector. 
%
% Optional key value pairs:
%   None.
%
% See also:
%   ObserverVecToParams, ComputeCIEConeFundamentals
%

% History:
%   08/09/19  dhb  Wrote it, because I have to do one fun thing this summer.

x = zeros(8,1);
x(1) =  params.indDiffParams.dlens;
x(2) =  params.indDiffParams.dmac;
x(3) =  params.indDiffParams.dphotopigment(1);
x(4) =  params.indDiffParams.dphotopigment(2);
x(5) =  params.indDiffParams.dphotopigment(3);
x(6) =  params.indDiffParams.lambdaMaxShift(1);
x(7) =  params.indDiffParams.lambdaMaxShift(2);
x(8) =  params.indDiffParams.lambdaMaxShift(3);
