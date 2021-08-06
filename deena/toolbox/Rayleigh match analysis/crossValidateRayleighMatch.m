function meanCrossValErr = crossValidateRayleighMatch(primaryMatchSpds,refMatchSpds,...
    lowerBounds,upperBounds,AEq,BEq,nOverallRuns,varargin)
% Computes the cross-validated error when parameters are fit to observers'
% Rayleigh matching data using various models
%
% Syntax:
%   crossValidateRayleighMatch(primaryMatchSpds,refMatchSpds,...
%   lowerBounds,upperBounds,AEq,BEq,nOverallRuns)
%
% Description:
%    Computes the cross-validated error associated with various parameter
%    fits to observers' Rayleigh matching data on the OneLight. Takes in
%    arrays of (radiometer-measured) primary and reference spds which were
%    identified by the observer as matches, as well as arrays containing
%    parameter constraints for the various models of interest. Performs
%    cross-validation for each model, with the number of iterations
%    determined by the number of times the matches were repeated at a given
%    reference wavelength during the experiment (in each iteration, a
%    different match is left out for each reference wavelength). Repeats
%    this procedure for the number of runs specified by nOverallRuns, and
%    averages the cross-validated fit error for each model across all runs.
%
% Inputs:
%    primaryMatchSpds - 1 x nReferenceLights cell array of (radiometer 
%                       measured) primary mixture spds which an observer  
%                       identified as Rayleigh matches. Each entry contains
%                       spectra for a different reference wavelength,
%                       entered as an spdLength x nMatches matrix.
%    refMatchSpds     -1 x nReferenceLights cell array of (radiometer 
%                      measured) primary mixture spds which an observer  
%                      identified as Rayleigh matches. Each entry contains
%                      spectra for a different reference wavelength,
%                      entered as an spdLength x nMatches matrix.
%    lowerBounds      -nModels x 8 matrix. Each row contains the lower
%                      parameter bounds for a different fitting model. See
%                      findObserverParameters for details on parameters.
%    upperBounds      -nModels x 8 matrix. Each row contains the upper
%                      parameter bounds for a different fitting model. See
%                      findObserverParameters for details on parameters.
%    AEq              -Cell array with (nModels) elements. Each entry contains
%                      the linear equality constraints on parameters for a
%                      different fitting model ([] = no constraint).
%    BEq              -Cell array with (nModels elements). Each entry
%                      contains the linear equality sum for a different
%                      fitting model ([] = no constraint).
%    nOverallRuns     -Number of times to run the overall cross-validation
%                      procedure. The cross-validated fit errors from each
%                      loop are averaged.
%
% Outputs:
%    crossValError  -Row vector where each entry is the cross-validated fit
%                    error for a different model of interest
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
%   06/30/21  dce, dhb  Wrote it.
%   08/05/21  dce       Edited to allow different numbers of matches for
%                       each reference wavelength 

% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('initialConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Generate base observer
baseObserver = genRayleighObserver('age',p.Results.age,'fieldSize',p.Results.fieldSize,...
    'coneVec',p.Results.initialConeParams,'opponentParams',p.Results.opponentParams);

% Number of models we are evaluating
nModelsToEval = length(BEq);

% Number of reference wavelengths provided
nRefWls = size(refMatchSpds);

% Matrix for storing cross-validated fit error from each model in each run
crossValErr = zeros(nOverallRuns,nModelsToEval);

% We track how many times the match at each reference wl was repeated.
nRepeats = zeros(1,nRefWls);
for rr = 1:nRefWls
    nRepeats(rr) = size(refMatchSpds{rr},2);
end 

% Run cross-validation for the specified number of runs 
for kk = 1:nOverallRuns
    % Set up cross-validation indices for the run 
    % Each entry of the cell array contains the order in which the matches
    % for that reference wavelength should be left out of the analysis. 
    leaveOneOutOrder = cell(1,nRefWls);
    for rr = 1:nRefWls
        leaveOneOutOrder{rr} = Shuffle(1:nRepeats(rr));
    end
    
    % Matrix to store error for each model and each iteration in the given
    % run
    crossValErrRun = zeros(nModelsToEval,minRepeats);
    
    % Perform cross validation. The number of iterations equals the minimum
    % number of repeats made for any reference wavelenth.
    for cc = 1:min(nRepeats)
        % Collect up the matches to fit for this iteration
        matchesToFit = cell(1,nRefWls);
        matchToEval = zeros(1,nRefWls);
        pSpdsToFit = [];
        refSpdsToFit = [];
        pSpdsToEval = zeros(spdLength,nRefWls);
        refSpdsToEval = zeros(spdLength,nRefWls);
        for rr = 1:nRefWls
            % We fit nRepeats(rr)-1 matches for each reference wavelength.
            % This code defines the indices of the ones being fit.
            matchesToFit{rr} = setdiff(1:nRepeats(rr),leaveOneOutOrder{rr}(cc));
            
            % We evaluate fit to one match for each reference wavelength.
            % This code defines the index of the one being fit.
            matchToEval(rr) = leaveOneOutOrder{rr}(cc);
            
            % Select the matches to fit and evaluate for the given reference
            % wavelength on this round of cross-validation
            pSpdsToFit = [pSpdsToFit,primaryMatchSpds{rr}(:,matchesToFit{rr})];
            refSpdsToFit = [refSpdsToFit,refMatchSpds{rr}(:,matchesToFit{rr})];
            
            pSpdsToEval(:,rr) = primaryMatchSpds{rr}(:,matchToEval(rr));
            refSpdsToEval(:,rr) = refMatchSpds{rr}(:,matchToEval(rr));
        end
        
        % Fit parameters using the specified matches and the provided bounds
        % for each of the models we are comparing, and evaluate the error for
        % the left-out matches.
        for mm = 1:nModelsToEval
            % Fit parameters using the selected matches
            [fitParams,~,~] = findObserverParameters(refSpdsToFit,pSpdsToFit,...
                'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
                'opponentParams',p.Results.opponentParams,'initialConeParams',...
                p.Results.initialConeParams,'minimizeConeErr',false,...
                'lowerBounds',lowerBounds(mm,:),'upperBounds',upperBounds(mm,:),...
                'AEq',AEq{mm},'BEq',BEq{mm},'errScalar',p.Results.errScalar);
            
            % Evaluate fit on the left-out data
            crossValErrRun(mm,cc) = findMatchError(fitParams,baseObserver,...
                refSpdsToEval,pSpdsToEval,'errScalar',1,'findConeErr',false);
        end
    end
    % For each model, average the fit error across cross validation 
    % iterations (sum of squares)
    for mm = 1:nModelsToEval
        crossValErr(kk,mm) = norm(crossValErrRun(mm,:));
    end
end
% Average cross-validated error from each model across all runs  
meanCrossValErr = mean(crossValErr,1);
end