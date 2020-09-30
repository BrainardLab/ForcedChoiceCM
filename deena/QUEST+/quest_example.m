% My attempt to code the Quest procedure for a toy example - finding the 
% mean of a cumulative normal distribution. This follows the 
% mQuestPlus notes from dhb 

% Define initial values for search 
cVec = linspace(0,1,100);  % Stimulus value can range from 0 to 1
muVec = linspace(0,1,100); % Set threshold so it can take on any possible stimulus value
sigma = 0.1; 

% Define the true distribution 
true_mu = 0.2; 
true_pdf = normcdf(cVec,true_mu,sigma);

% Define the prior 
p_mu = unifpdf(muVec,muVec(1),muVec(end))*1/100; % Uniform prior for mu
posterior = p_mu; % Posterior is initially set to the prior

% Arrays for storing data. 
posteriorsByOutcome = zeros(length(cVec),length(muVec),2); 
expected_entropies = zeros(length(cVec),1); 
nTrials = 64;              % Number of trials 
trials = zeros(nTrials,2); % Trial data with format [c, outcome]

% Compute outcome probabilities for each c-mu combination: p (r | muj, ci). 
% Arrays follow the convention [false, true]
outcomeProbs = zeros(length(cVec),length(muVec),2); 
outcomeProbsByC = zeros(length(muVec),2);  
for i = 1:length(cVec)
    outcomeProbs(i,:,1) = normcdf(cVec(i),muVec,sigma);
    outcomeProbs(i,:,2) = 1- outcomeProbs(i,:,1);
    % Is this calculation correct?
    outcomeProbsByC(i,:) = sum(squeeze(outcomeProbs(i,:,:).*p_mu))';  
end

for kk = 1:nTrials
    % Compute posteriors by outcome for each c-mu combination: p(mu j | ci, r). 
    % Each column represents the posteriors for a particular mu over various 
    % values of c for a given outcome. In each trial, the previous
    % posterior is used as the prior. 
    for i = 1:length(cVec)
        for j = 1:length(muVec)
            posteriorsByOutcome(i,j,:) = (outcomeProbs(i,j,:)....
                .*posterior(j))./sum(outcomeProbs(i,:,:).*posterior);
        end
    end
    
    % Compute expected entropies for each stimulus, then select c that has
    % minimum expected entropy.
    for i = 1:length(cVec)
        e_false = -nansum(posteriorsByOutcome(i,:,1).*log2(posteriorsByOutcome(i,:,1)));
        e_true = -nansum(posteriorsByOutcome(i,:,2).*log2(posteriorsByOutcome(i,:,2)));
        expected_entropies(i) = e_false*outcomeProbsByC(i,1)+e_true*outcomeProbsByC(i,2);
        [~,minC] = min(expected_entropies);
    end
    
    % Sample outcome for the trial
    trialOutcomeProportion = [normcdf(cVec(minC),true_mu,sigma),...
        1- normcdf(cVec(minC),true_mu,sigma)];
    outcomeVector = mnrnd(1,trialOutcomeProportion);
    outcome = find(outcomeVector);
    
    % Update posterior given the outcome
    trials(kk,:) = [minC,outcome];
    posterior = posteriorsByOutcome(minC,:,outcome);
end

% Find estimated threshold
[~,threshInd] = max(posterior); 
estThresh = muVec(threshInd); 

% Plot results
figure; 
hold on;
plot(cVec,true_pdf); 
plot(cVec(trials(:,1)), trials(:,2)-1, 'o ');
plot(true_mu,0.5,'m*');
plot(estThresh,0.5,'g*');
legend('True pdf','Trials','True Mean','Estimated Mean'); 
title('Quest Trial Sequence');
xlabel('Stimulus Value'); 
ylabel('Probability'); 