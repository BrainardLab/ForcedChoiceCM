function plotSampledRayleighMatch(nObservers,baseParams,paramsToVary,...
    p1,p2,test,method,varargin)
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
%    Produces a series of plots
%
% Optional key-value pairs:
%    'foveal'            -logical indicating whether we are making foveal
%                         matches. Default is true.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.01.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'monochromatic'     -When using predicted matches, logical indicating
%                         to use monochromatic spds insted of those
%                         generated by the OneLight. Default is false.
%    'nReversals'        -When using a simulated observer, number of
%                         reversals the observer must make before
%                         changing step size. Enter as a 2-element vector -
%                         the first element is the number of reversals for
%                         intermediate step sizes, the second is the number
%                         needed for the smallest step size. Default is
%                         [1 4].
%    'nBelowThreshold'   -When using a simulated observer with
%                         threshold matching, number of pairs below
%                         threshold required before recording a match.
%                         Default is 1.
%    'thresholdScaleFactor' -When using a simulated observer with
%                            threshold matching, scale factor for matching
%                            threshold. Default is 0.5.

% History:
%   07/06/20  dce       Wrote it.

% Parse input
p = inputParser;
p.addParameter('foveal',true,@(x)(islogical(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.01,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('monochromatic',false,@(x)(islogical(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Standard deviation of parameters (Asano 2015)
sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];

% Data-storing arrays. Each row represents a different observer's params
sampledParams = [];
recoveredParams = [];
testIntensitiesSim = [];
testIntensitiesPred = [];
primaryRatiosSim = [];
primaryRatiosPred = [];

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
    [testSpds,primarySpds,testIntensitiesSimObs,primaryRatiosSimObs] = ...
        getMatchSeries(subjID,observerParams,p1,p2,test,method,'foveal',...
        p.Results.foveal,'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'nObserverMatches',...
        p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,'thresholdScaleFactor',...
        p.Results.thresholdScaleFactor,'rayleighPlots',false);
    testIntensitiesSim = [testIntensitiesSim;testIntensitiesSimObs];
    primaryRatiosSim = [primaryRatiosSim;primaryRatiosSimObs];
    
    % Recover observer parameters
    if p.Results.foveal
        fieldSize = 2;
    else
        fieldSize = 10;
    end
    [calcParams,error] = findObserverParameters(testSpds,primarySpds,...
        'age',p.Results.age,'fieldSize',fieldSize);
    recoveredParams = [recoveredParams;calcParams];
    
    % Recover predicted test intensity and primary ratio based on
    % recovered parameters
    [~,~,testIntensitiesPredObs,primaryRatiosPredObs] = ...
        getMatchSeries([subjID,'_pred'],[calcParams 0],p1,p2,test,...
        'predicted','foveal',...
        p.Results.foveal,'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic);
    testIntensitiesPred = [testIntensitiesPred;testIntensitiesPredObs];
    primaryRatiosPred = [primaryRatiosPred;primaryRatiosPredObs];
end

%% Plots
coneParamNames = {'Lens Density','Macular Pigment Density',...
    'L Photopigment Density','M Photopigment Density',...
    'S photopigment density','L Lambda Max','M Lambda Max', 'S Lambda Max'};

% Actual vs predicted plots for each parameter
for k = 1:length(paramsToVary)
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
    legend('Parameters','y=x');
end

%% Cone spectral sensitivities
for k = 1:nObservers
    S = [380 2 201];
    sampledObserver = genRayleighObserver('coneVec',...
        [sampledParams(k,:) 0],'S',S);
    recoveredObserver = genRayleighObserver('coneVec',...
        [recoveredParams(k,:) 0],'S',S);
    
    figure();
    plot(SToWls(S),sampledObserver.T_cones(1:2,:),'b-',SToWls(S),...
        recoveredObserver.T_cones(1:2,:),'r-');
    legend('sampled observer','recovered observer'); % need to fix legend
    theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g', k);
    title(theTitle);
end

%% Generalized Pitt diagram
for i = 1:nObservers
    figure();
    plot(primaryRatiosSim(i,:),testIntensitiesSim(i,:),'b-o',...
        primaryRatiosPred(i,:),testIntensitiesPred(i,:),'r-o');
    theTitle = sprintf('Generalized Pitt Diagram, observer %g',i);
    title(theTitle);
    xlabel('Primary Ratio');
    ylabel('Test Intensity');
    legend('Simulated','Predicted');
end

end