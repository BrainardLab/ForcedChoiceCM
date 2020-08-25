function sampledParams = sampleRayleighObservers(nObservers,baseParams,...
    paramsToVary)
% Samples a set of observers for Rayleigh matching 
%
% Syntax:
%   sampleRayleighIbservers(nObservers,baseParams,paramsToVary);
%
% Description:
%    Takes in a vector of cone individual difference parameters as a
%    starting point, as well as information about which parameters to vary.
%    Creates a series of n simulated observers by sampling the chosen
%    parameters around their standard deviations.
%
% Inputs:
%    nObservers     -Number of simulated observers to sample.
%    baseParams     -Eight-element numeric vector of individual difference
%                    parameters (see ObserverVecToParams for description).
%                    The values are used as means for observer sampling.
%    paramsToVary   -Eight-element numeric vector of ones and zeros
%                    indicating which individual difference parameters
%                    should be varied (the noise parameter is excluded).
%                    Parameters set to 1 will be sampled around their
%                    standard deviation, while parameters set to 0 will
%                    stay at the values specified in baseParams.
% Outputs:
%    sampledObservers  -nObservers x 8 array of sampled observer parameters.
%
% Optional key-value pairs:
%    none

% History 
%   07/29/20  dce      Separated from sampleRayleighMatch.

% Standard deviation of parameters (Asano 2015). 1:5 are percent 
% deviations, while 6:8 are deviations in nm.
sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];

% Output array. Each row is a different observer
sampledParams = zeros(nObservers,length(baseParams));

% Sample parameters for each observer and add to array 
for i = 1:nObservers
    observerParams = baseParams;
    for j = 1:length(sds)
        if logical(paramsToVary(j))
            observerParams(j) = normrnd(baseParams(j),sds(j));
        end
    end
    sampledParams(i,:) = observerParams;
end 