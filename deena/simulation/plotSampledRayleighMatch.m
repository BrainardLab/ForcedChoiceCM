function plotSampledRayleighMatch(nObservers,baseParams,paramsToVary,p1,p2,test,method)
% Tests the OLRayleighMatch simulation and makes associated sampling plots. 

% Syntax:
%   plotSampledRayleighMatch(nObservers,baseParams,paramsToVary,p1,p2,test,method)
%
% Description:
%    Takes in a vector of cone individual difference parameters as a 
%    starting point, as well as information about which parameters to vary.
%    Creates a series of simulated observers by sampling the chosen 
%    parameters based on their standard deviations. For each observer, runs 
%    a Rayleigh matching simulation with a set of primary/test lights 
%    (OLRayleighMatch or computePredictedRayleighMatch) and uses the match
%    results to recover the individual difference parameters 
%    (findObserverParameters). 
%
%    The program produces three types of plots to validate how similar the
%    recovered parameters are to the sampled (true) parameters. The first 
%    is a plot of sampled vs predicted values for each variable parameter,
%    aggregated across all observers. The second is a comparative plot of 
%    cone spectral sensitivities based on the sampled and recovered 
%    parameters for each observer. The third, a "generalized Pitt diagram,"
%    plots primary ratio and test intensity settings that were found in the 
%    simulation for each observer, comparing them to the settings predicted 
%    from the recovered parameters. 
%
%    For both the second and third plot types, one figure is produced for
%    each observer. For this reason, these two plot types are optional, and
%    it is generally best not to include them when sampling a large number
%    of observers. 
%
% Inputs:
%    nObservers     -Number of simulated observers to sample.
%    baseParams     -Nine-element numeric vector of individual difference 
%                    parameters (see ObserverVecToParams for description). 
%                    The values are used as means for observer sampling.
%    paramsToVary   -Eight-element numeric vector of ones and zeros
%                    indicating which individual difference parameters 
%                    should be varied (the noise parameter is excluded). 
%                    Parameters set to 1 will be sampled around their 
%                    standard deviation, while parameters set to 0 will 
%                    stay at the values specified in baseParams.
%    p1             -Integer or numeric vector of desired wavelengths
%                    for the first primary.
%    p2             -Integer or numeric vector of desired wavelengths
%                    for the second primary.
%    test           -Integer or numeric vector of desired wavelengths
%                    for the test light.
%    method         -Character vector for match method. Choose either
%                    'predicted', 'forcedChoice', or 'threshold'.
%
% Outputs:
%    none
%
% Optional key-value pairs: 

% History:
%   07/06/20  dce       Wrote it.

% Standard deviation of parameters (Asano 2015)
sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];

% Data-storing arrays. Each row represents a different observer's params
sampledParams = [];
recoveredParams = [];

% For each observer: sample parameters, make matches, and use matches to
% recover parameters
for i = 1:nObservers
    % Sample parameters
    observerParams = baseParams;
    for j = 1:length(sds)
        if paramsToVary(j) ~= 0
            observerParams(j) = normrnd(baseParams(j),sds(j));
        end
    end
    sampledParams = [sampledParams;observerParams];
    subjID = mat2str(observerParams);
    
    % Make a series of Rayleigh matches for the observer
    [testSpds,primarySpds] = ...
        getMatchSeries(subjID,observerParams,p1,p2,test,method);
    
    % Recover observer parameters
    [calcParams,error] = findObserverParameters(testSpds,primarySpds);
    recoveredParams = [recoveredParams;calcParams];
end

%% Plots
coneParamNames = {'Lens Density','Macular Pigment Density',...
    'L Photopigment Density','M Photopigment Density',...
    'S photopigment density','L Lambda Max','M Lambda Max', 'S Lambda Max'};

% Actual vs predicted plots for each varied parameter
for k = 1:length(paramsToVary)
    if paramsToVary(k) ~= 0
        xVals = sampledParams(:,k);   % Predicted parameters
        yVals = recoveredParams(:,k); % Recovered params
        figure();
        hold on;
        plot(xVals,yVals,'b* ');
        refline(1,0);
        theTitle = sprintf('%s Predicted vs Actual',cell2mat(coneParamNames(k)));
        title(theTitle);
        xlabel('Sampled Parameters');
        ylabel('Recovered Parameters');
        % Add legend
    end
end

% Cone spectral sensitivities
for k = 1:nObservers
    % figure out what S is
    wls = SToWls(S);
    sampledObserver = genRayleighObserver('coneVec',...
        [sampledParams(k,:) 0],'S',S);
    recoveredObserver = genRayleighObserver('coneVec',...
        [recoveredParams(k,:) 0],'S',S);
    figure();
    plot(wls,sampledObserver.T_cones,'b-',wls,recoveredObserver.T_cones,'r-');
    legend('sampled observer','recovered observer');
    theTitle = sprintf('Cone Spectral Sensitivities, observer %g', k);
    title(theTitle);
end

% Generalized Pitt diagram
nMatches = length(p1)*length(p2)*length(test);
for i = 1:nObservers
    % Figure out how to recover lambdas
    figure(); 
    lambdaSims = []; 
    lambdaPreds = []; 
    testIntensitySims = []; 
    testIntensityPreds = [];
    for j = 1:nMatches
        [~,~,lambdaPred,testIntensityPred] =computePredictedRayleighMatch...
            (p1,p2,test,[recoveredParams(i,:),0]);
        lambdaPreds = [lambdaPreds,lambdaPred]; 
        testIntensityPreds = [testIntensityPreds,testIntensityPred]; 
    end
    plot(lambdaSims,testIntensitySims,'b-o',lambdaPreds,...
        testIntensityPreds,'r-o'); 
    theTitle = sprintf('Generalized Pitt Diagram, observer %g',i);
    title(theTitle);
    xlabel('Primary Ratio'); 
    ylabel('Test Intensity'); 
    legend('Simulated','Predicted'); 
end

end