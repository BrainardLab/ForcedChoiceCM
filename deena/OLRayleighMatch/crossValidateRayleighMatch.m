function crossValErr = crossValidateRayleighMatch(primaryMatchSpds,refMatchSpds,...
    lowerBounds,upperBounds,AEq,BEq,varargin)
% Computes the cross-validated error associated with a set of Rayleigh 
% matches with parameters fit using different models
%
% Syntax:
%   crossValidateRayleighMatch(primaryMatchSpds,refMatchSpds,...
%   lowerBounds,upperBounds,AEq,BEq)
%
% Description:
%    Given a set of primary and reference spds which an observer identified
%    as Rayleigh matches, performs cross-validation to determine which of
%    several models best accounts for the observer's data. 
%
% Inputs:
%    primaryMatchSpds -[spdLength x nCrossValIters x nReferenceLights]
%                      array of (radiometer-measured) primary mixture spds 
%                      which an observer identified as Rayleigh matches
%    refMatchSpds     -[spdLength x nCrossValIters x nReferenceLights]
%                      array of (radiometer-measured) reference spds 
%                      which an observer identified as Rayleigh matches
%    lowerBounds      -nModels x 8 matrix. Each row contains the lower
%                      parameter bounds for a different fitting model. See
%                      findObserverParameters for details on parameters.
%    upperBounds      -nModels x 8 matrix. Each row contains the upper
%                      parameter bounds for a different fitting model. See
%                      findObserverParameters for details on parameters.
%    AEq              -nModels x 8 matrix. Each row contains the linear 
%                      equality constraints on parameters for a different 
%                      fitting model.
%    BEq              -Vector with (nModels) elements containing the linear 
%                      equality sum for each model tested.
%                                         
% Outputs:
%    crossValError  -Row vector where each entry is the cross-validated fit
%                    error for a different model of interest
%
% Optional key-value pairs:
%    'age'           -Observer age in years. Default is 32.
%    'fieldSize'     -Field size in degrees. Default is 2.
%    'initialConeParams' -1x8 numerical vector of cone individual  
%                         difference parameters. Default is zeros(1,8);
%    'opponentParams'    -1x4 numerical vector of opponent contrast
%                         parameters. Default is [40.3908 205.7353 62.9590 1.0000].
%    'errScalar'         -Integer for scaling the match error in the 
%                         fitting program, to improve search. Default is 100.

% History:
%   06/30/21  dce, dhb  Wrote it.

% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('initialConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.Parse(varargin{:});

% Generate base observer
baseObs = genRayleighObserver('age',p.Results.age,'fieldSize',p.Results.fieldSize,...
    'coneVec',p.Results.initialConeParams,'opponentParams',p.Results.opponentParams);

% Set up cross-validation indices
% Each row of the matrix set up here is the order to leave one
% measurement out, randomized separately for each reference light.
[spdLength,nCrossValIters,nReferenceLights] = size(refMatchSpds);
leaveOneOutOrder = zeros(nReferenceLights,nCrossValIters);
for rr = 1:nReferenceLights
    leaveOneOutOrder(rr,:) = Shuffle(1:nCrossValIters);
end

nModelsToEval = length(BEq);  % Number of models we are evaluating
crossValError = zeros(nModelsToEval,nCrossValIters); % Store error for each model 

% Loop through the reference lights, and perform cross validation for the
% specified number of iterations
for cc = 1:nCrossValIters   
    % Collect up the matches to fit for this iteration
    matchesToFit = zeros(nReferenceLights,nCrossValIters-1);
    matchToEval = zeros(nReferenceLights,1);
    pSpdsToFit = [];
    refSpdsToFit = [];
    pSpdsToEval = zeros(spdLength,nRefWls);
    refSpdsToEval = zeros(spdLength,nRefWls);
    for rr = 1:nReferenceLights
        % We fit nCrossValIters-1 matches for each reference wavelength.
        % This code defines the indices of the ones being fit.
        matchesToFit(rr,:) = setdiff(1:nCrossValIters,leaveOneOutOrder(rr,cc));
        
        % We evaluate fit to one match for each reference wavelength.
        % This code defines the index of the one being fit.
        matchToEval(rr) = leaveOneOutOrder(rr,cc);
        
        % Select the matches to fit and evaluate for the given reference
        % wavelength on this round of cross-validation
        pSpdsToFit = [pSpdsToFit,primaryMatchSpds(:,matchesToFit(rr,:),rr)];
        refSpdsToFit = [refSpdsToFit,refMatchSpds(:,matchesToFit(rr,:),rr)];
        
        pSpdsToEval(:,rr) = primaryMatchSpds(:,matchToEval(rr),rr);
        refSpdsToEval(:,rr) = refMatchSpds(:,matcheToEva(rr),rr);
    end
    
    % Fit parameters using the specified matches and the provided bounds
    % for each of the models we are comparing, and evaluate the error for 
    % the remaining match. 
    for mm = 1:nModelsToEval
        % Fit parameters using the selected matches
        [fitParams,~,~] = findObserverParameters(refSpdsToFit,pSpdsToFit,...
            'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
            'opponentParams',p.Results.opponentParams,'initialConeParams',...
            p.Results.initialConeParams,',minimizeConeErr',false,...
            'lowerBounds',lowerBounds(mm,:),'upperBounds',upperBounds(mm,:),...
            'AEq',AEq(mm,:),'BEq',BEq(mm),'errScalar',p.Results.errScalar,...
            'minimizeConeErr',false);
        
        % Evaluate fit on the left-out data
        crossValError(mm,cc) = findObserverError(fitParams,baseObserver,...
            refSpdsToEval,pSpdsToEval,'errScalar',1,'findConeErr',false);
    end
end
% For each model, average the fit error across cross validation iterations
crossValErr = zeros(1,nModelsToEval);
for mm = 1:nModelsToEval
    crossValErr(mm) = sqrt(sum(crossValError(mm,:).^2));
end 
end 