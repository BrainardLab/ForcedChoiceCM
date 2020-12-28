function paramErr = findParamRecoveryError(observerParams1,observerParams2)
% Find root mean square error between sets of cone parameter
%
% Syntax:
%   findParamRecoveryError(observerParams1,observerParams2)
%
% Description:
%    Takes two sets of observer parameters and computes the root mean
%    square error of the differences between the parameters, indicating how
%    close the parameter sets are to one another.
%
% Inputs:
%    observerParams1  -Vector of eight individual difference parameters for 
%                      the first observer. See ObserverVecToParams for a 
%                      full description.
%    observerParams1  -Vector of eight individual difference parameters for 
%                      the second observer. See ObserverVecToParams for a 
%                      full description.
%
% Outputs:
%    paramErr    -Root mean square error of the parameters
%
% Optional key-value pairs:
%    none
%                        
% History 
%   12/28/20  dce       Wrote it. 

if length(observerParams1)~=length(observerParams2)
    error('Input parameter vectors have unequal lengths');
end
squaredDiffs = (observerParams1-observerParams2).^2;
paramErr = sqrt(mean(squaredDiffs));
end