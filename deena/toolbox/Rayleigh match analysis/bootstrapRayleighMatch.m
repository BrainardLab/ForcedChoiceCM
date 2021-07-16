function [confidenceIntervals,bootstrapFitParams] = ...
    bootstrapRayleighMatch(primaryMatchSpds,refMatchSpds,...
    lowerBounds,upperBounds,AEq,BEq,nIterations,varargin)
% Bootstraps the parameter fit of a specified model to observers' Rayleigh
% matching data
%
% Syntax:
%   bootstrapRayleighMatch(primaryMatchSpds,refMatchSpds,...
%   lowerBounds,upperBounds,AEq,BEq,nIterations)
%
% Description:
%    Uses bootstrapping to find confidence intervals associated with fitting
%    a particular version of the Asano model to observers' Rayleigh matching
%    data. Takes in arrays of (radiometer-measured) primary and reference 
%    spds which were identified by the observer as matches, as well as 
%    arrays containing parameter constraints for the model of interest. 
%    Performs boostrapping for a specified number of iterations- on each 
%    iteration, the program resamples the primary and reference match spds
%    and uses the resampled spectra to fit the model of interest. Computes
%    confidence intervals based on the bootstrapped parameter fits. 
%
% Inputs:
%    primaryMatchSpds -[spdLength x nMatchesPerRef x nReferenceLights]
%                      array of (radiometer-measured) primary mixture spds
%                      which an observer identified as Rayleigh matches
%    refMatchSpds     -[spdLength x nMatchesPerRef x nReferenceLights]
%                      array of (radiometer-measured) reference spds
%                      which an observer identified as Rayleigh matches
%    lowerBounds      -8-element vector with lower parameter limits for 
%                      model fit. See findObserverParameters for details on
%                      parameters.
%    upperBounds      -8-element vector with upper parameter limits for 
%                      model fit. See findObserverParameters for details on
%                      parameters.
%    AEq              -vector containing linear equality constraints on 
%                      parameters for the model fit ([] = no constraint).
%    BEq              -Vector containing the linear equality sum for the
%                      model fit ([] = no constraint).
%    nIterations      -Number of times to run the bootstrapping procedure.
%
% Outputs:
%
% Optional key-value pairs:
%    'age'               -Observer age in years. Default is 32.
%    'fieldSize'         -Field size in degrees. Default is 2.
%    'initialConeParams' -1x8 numerical vector of cone individual
%                         difference parameters. Default is zeros(1,8);
%    'opponentParams'    -1x4 numerical vector of opponent contrast
%                         parameters. Default is [40.3908 205.7353 62.9590 1.0000].
%    'errScalar'         -Integer for scaling the match error in the
%                         fitting program, to improve search. Default is 100.

% History:
%   07/02/21  dce      Wrote it.
%   07/09/21  dce      Added confidence interval computation

% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('initialConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Check the size of our input array. The first dimension is the spd length,
% the second is the number of matches for each reference wavelength, and
% the third is the number of reference wavelengths tested
[~,nMatchesPerRef,nRefWls] = size(refMatchSpds);

% Matrix for storing parameters fit from each iteration
nParams = 8; % The Asano model has 8 parameters
bootstrapFitParams = zeros(nIterations,nParams);

% Bootstrapping for the specified number of iterations 
for kk = 1:nIterations
    % For each reference wavelength, we sample nMatchesPerRef matches for
    % analysis in this iteration, with replacement. Here, we determine  
    % which indices will be sampled for each wavelength in this iteration.
    trialInds = randi(nMatchesPerRef,nRefWls,nMatchesPerRef);
    
    % Collect the spds which will be used for fitting here
    trialRefSpds = [];
    trialPrimarySpds = [];
    for rr = 1:nRefWls
        for mm = 1:nMatchesPerRef
            trialRefSpds = [trialRefSpds,refMatchSpds(:,trialInds(rr,mm),rr)];
            trialPrimarySpds = [trialPrimarySpds,primaryMatchSpds(:,trialInds(rr,mm),rr)];
        end 
    end 
   
    % Fit parameters using the specified matches and the provided bounds 
    [bootstrapFitParams(kk,:),~,~] = ...
        findObserverParameters(trialRefSpds,trialPrimarySpds,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'opponentParams',p.Results.opponentParams,'initialConeParams',...
        p.Results.initialConeParams,'minimizeConeErr',false,...
        'lowerBounds',lowerBounds,'upperBounds',upperBounds,...
        'AEq',AEq,'BEq',BEq,'errScalar',p.Results.errScalar);
end
    
% Compute confidence intervals based on fit parameters - central 65% of
% bootstrapped distribution. Lower bounds are in the first row, and upper
% bounds are in the second row.
confidenceIntervals = zeros(2,8);
confidenceIntervals(1,:) = prctile(bootstrapFitParams,17.5,1);
confidenceIntervals(2,:) = prctile(bootstrapFitParams,82.5,1);
end