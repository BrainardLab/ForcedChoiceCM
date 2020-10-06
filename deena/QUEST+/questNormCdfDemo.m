% questNormCdfDemo  
% Demonstrates basic use of the QUEST procedure
%
% Description:
%    This script illustrates the Quest procedure using a simple example - 
%    finding the mean of a cumulative Normal distribution. All core 
%    functionality is included in the script itself rather than external 
%    functions. The code is intended for demonstration purposes, and it is  
%    not optimized for speed. It follows the mQuestPlus notes from dhb. 

% History: 
%    09/29/20  dce  - Wrote it, based on the mQUESTPlus code

% Define initial values 
cVec = linspace(0,1,100);  % Possible values for stimulus intensity, ranging from 0 to 1 
muVec = linspace(0,1,100); % Possible parameter values for the mean
sigma = 0.1;               % Known standard deviation

% Define the actual mean and the resulting cumulative Normal distribution 
true_mu = 0.2; 
true_pdf = normcdf(cVec,true_mu,sigma); 

% Define the prior. The first trial uses the uniform prior p_mu, but
% this is updated for later trials.
p_mu = unifpdf(muVec,muVec(1),muVec(end))*1/100; % Uniform prior for mu
prior = p_mu; 

% % Compute the probabilities by outcome for each c-mu combination:  
% % p (r |c_i,mu_j). The third dimension of arrays follow the convention 
% % [P(false), P(true)]. 
outcomeProbs = zeros(length(cVec),length(muVec),2); 
outcomeProbsByC = zeros(length(muVec),2);  
for i = 1:length(cVec)
    % Calculate probabilities for a given value of c at each value of mu
    outcomeProbs(i,:,1) = 1 - normcdf(cVec(i),muVec,sigma); % P(false)
    outcomeProbs(i,:,2) = normcdf(cVec(i),muVec,sigma);     % P(true)
    
    % Marginalize the outcome probabilities over mu to get the outcome
    % probabilities for a given stimulus value c_i. 
    outcomeProbsByC(i,:) = sum(squeeze(outcomeProbs(i,:,:).*p_mu))';
end

% Initialize arrays for storing data. 
nTrials = 64;              % Number of trials 
trials = zeros(nTrials,2); % Trial data with format [c, outcome]
posteriorsByOutcome = zeros(length(cVec),length(muVec),2); 
expected_entropies = zeros(length(cVec),1); 
% 
% % Set up Quest+ functions for testing
% questData = qpInitialize('stimParamsDomainList',{0:0.0101:1}, ...
%     'psiParamsDomainList',{0:0.0101:1, 0.1, 0},'qpPF',@qpPFNormal);
% simulatedObserverFun = @(x) qpSimulatedObserver(x,@qpPFNormal,[true_mu, 0.1 0]);

% Loop through trials 
for kk = 1:nTrials
    % Compute posteriors by outcome for each c-mu combination: p(mu_j|c_i,r). 
    % Values of c are in the row index, values of mu are in the column 
    % index, and the third dimension holds outcome probabilities 
    % [P(false),P(true)]. 
    for i = 1:length(cVec)
        for j = 1:length(muVec)
            posteriorsByOutcome(i,j,:) = (outcomeProbs(i,j,:)....
                .*prior(j))./sum(outcomeProbs(i,:,:).*prior);
        end

    end
    
    % Compute expected entropies for each stimulus. 
    % e_false is the entropy for a false response, e_true is the entropy
    % for a true response, and expected_entropy is an average of the two
    % weighted by their probabilities.
    for i = 1:length(cVec)
        e_false = -nansum(posteriorsByOutcome(i,:,1).*log2(posteriorsByOutcome(i,:,1)));
        e_true = -nansum(posteriorsByOutcome(i,:,2).*log2(posteriorsByOutcome(i,:,2)));
        expected_entropies(i) = e_false*outcomeProbsByC(i,1)+e_true*outcomeProbsByC(i,2); 
    end
    
    % Select the value of c which leads to the minimum expected entropy. 
    % This will be the stimulus value for the next trial. 
    [~,minCInd] = min(expected_entropies);
    
    % Using the chosen value of c, simulate a trial by sampling from a 
    % multinomial distribution with outcome proportions [P(false), P(true)] 
    trialOutcomeProportion = [1 - normcdf(cVec(minCInd),true_mu,sigma),...
        normcdf(cVec(minCInd),true_mu,sigma)];  % [P(false), P(true)]
    outcomeVector = mnrnd(1,trialOutcomeProportion); 
    outcome = find(outcomeVector);               
    
    % Now that the outcome has been determined, we know which of the
    % posteriors by outcome represents the actual posterior. This posterior
    % is used as the prior in subsequent trials. 
    posterior = posteriorsByOutcome(minCInd,:,outcome);
    prior = posterior; 
    trials(kk,:) = [cVec(minCInd),outcome]; % Update trials data 
    
%     % For testing
%     qStim = qpQuery(questData);
%     qOutcome = simulatedObserverFun(qStim);
%     questData = qpUpdate(questData,qStim,qOutcome);
%     qPrior = questData.posterior; 
end

% Once all trials are completed, estimate the threshold as the maximum of
% the posterior. 
[~,threshInd] = max(posterior); 
estThresh = muVec(threshInd); 

% Plot results. Trials are shown as open circles, with x positions
% representing the value of c. Their y positions are either 0 (false) or 1
% (true). 
figure; 
hold on;
plot(cVec,true_pdf); 
plot(trials(:,1), trials(:,2)-1, 'o ');
plot(true_mu,0.5,'m*');
plot(estThresh,0.5,'g*');
legend('True pdf','Trials','True Mean','Estimated Mean'); 
title('Quest Trial Sequence');
xlabel('Stimulus Value'); 
ylabel('Probability'); 