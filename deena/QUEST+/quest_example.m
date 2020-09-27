% My attempt to code the Quest procedure for a toy example - finding the 
% mean of a cumulative normal distribution. This follows the 
% mQuestPlus notes from dhb 

% Define initial values for search 
cVec = linspace(0,1,100);  % Stimulus value can range from 0 to 1
muVec = linspace(0,1,100); % Set threshold so it can take on any possible stimulus value
sigma = 0.1; 
nTrials = 64; % number of trials 

% Define the true distribution 
true_mu = 0.3; 
true_pdf = normcdf(cVec,true_mu, sigma);

% Define initial arrays. The third dimension follows the convention (false, true_ 
p_mu = unifpdf(muVec,muVec(1),muVec(end)); % Uniform prior for threshold
outcomeProbs = zeros(length(muVec),length(cVec),2);  
likelihoods = zeros(length(muVec),length(cVec),2);
posteriors = zeros(length(muVec),length(cVec),2); 
entropies = zeros(length(muVec),length(cVec),2); 

% Compute outcome probabilities for each c-mu combination 
for i = 1:length(muVec)
    for j = 1:length(cVec)
        outcomeProbs(i,j,:) = [normcdf(cVec(j),muVec(i),sigma),1- normcdf(cVec(j),muVec(i),sigma)]; 
    end 
end
% Compute posteriors for each c-mu combination
for i = 1:length(muVec)
    for j = 1:length(cVec)
        posteriors(i,j,:) = (outcomeProbs(i,j,:).*p_mu(i))./sum(outcomeProbs(:,j,:).*p_mu');
    end 
end
% Compute expected posterior entropy for each c-mu combination 
for i = 1:length(muVec)
    for j = 1:length(cVec)
        entropies(i,j,:) = 
    end 
end

for i = 1:nTrials
    % which trial should we be looking at?
    
    % Sample outcome for the trial
    
    % Update arrays
end 

