function [sampledParams,recoveredParams,error,rayleighSettings] = ...
    sampleRayleighMatch(nObservers,baseParams,paramsToVary,...
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
%    (findObserverParameters). Then, conducts a few basic analyses to
%    validate how similar the recovered parameters are to the sampled 
%    (true) parameters. 
%
%    The program produces four types of plots. The first is a plot of 
%    sampled vs predicted values for each parameter, aggregated across all 
%    observers. The second is a bar plot of the difference between the cone
%    sensitivities of the sampled and recovered parameters for each 
%    observer (root mean square error), compared to the difference for 
%    the sampled parameters and the standard observer. The third is a
%    comparative plot of cone spectral sensitivities based on the sampled 
%    and recovered parameters for each observer. The third, a "generalized 
%    Pitt diagram," plots primary ratio and test intensity settings that 
%    were found in the simulation for each observer, comparing them to the 
%    settings predicted from the recovered parameters.
%
%    By default, the third and fourth plot types are produced for the best
%    and the worst observers (identified based on how different their cone
%    spectral sensitivities are). However, there is also an option to 
%    produce the plots for every observer.
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
%    sampledParams    -nObservers x 9 array of the sampled individual
%                      difference parameters.
%    recoveredParams  -nObservers x 8 array of the recovered individual
%                      difference parameters.
%    error            -Struct with data on the cone spectral sensitivity
%                      error. Fields include coneErr (sampled vs recovered
%                      excitations), coneStandardErr (sampled vs standard
%                      observer excitations), coneErr540, and
%                      coneStandardErr540 (same as the previous two, but
%                      calculated for wavelengths over 540nm).
%    rayleighSettings -Struct with data on Rayleigh match settings. Fields
%                      include testIntensitiesSim (test intensities based
%                      on calculated parameters), testIntensitiesPred (test
%                      intensities predicted by recovered parameters),
%                      primaryRatiosSim (primary ratios based on calculated
%                      parameters), and primaryRatiosPred (primary ratios
%                      predicted by recovered parameters).
%
% Optional key-value pairs:
%    'makeAllPlots'      -logical indicating whether we are plotting the
%                         cone spectral sensitivities and generalized Pitt
%                         diagrams for each observer. Otherwise, plots are
%                         made for the best and the worst observers.
%                         Default is false.
%    'fieldSize          -Integer field size, in degrees. Default is 2.
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
%   07/10/20  dce       Edited plotting and added rms error output.

% Close stray figures
close all;

% Parse input
p = inputParser;
p.addParameter('makeAllPlots',false,@(x)(islogical(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
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

% Data-storing arrays. Each row holds a different observer's params.
sampledParams = [];           % Sampled cone parameters
recoveredParams = [];         % Cone parameters recovered by simulation
testIntensitiesSim = [];      % Calculated test intensities
testIntensitiesPred = [];     % Test intensities predicted from recovered params 
primaryRatiosSim = [];        % Calculated primary ratios
primaryRatiosPred = [];       % Primary ratios predicted from recovered params

% Arrays for spectral sensitivity error data, both across the gamut and at 
% wavelengths greather than 540nm. The first two hold rms error for
% the sampled parameters relative to the recovered parameters, while the
% second two hold rms error for the sampled parameters relative to the
% standard observer. 
coneErr = zeros(nObservers,1); 
coneErr540 = zeros(nObservers,1); 
coneStandardErr = zeros(nObservers,1);  
coneStandardErr540 = zeros(nObservers,1);  

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
    subjID = 'test_series';
    
    % Make a series of Rayleigh matches for the observer
    [testSpds,primarySpds,testIntensitiesSimObs,primaryRatiosSimObs] = ...
        getMatchSeries(subjID,observerParams,p1,p2,test,method,'fieldSize',...
        p.Results.fieldSize,'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'nObserverMatches',...
        p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,'thresholdScaleFactor',...
        p.Results.thresholdScaleFactor,'rayleighPlots',false,...
        'saveResults',false);
    testIntensitiesSim = [testIntensitiesSim;testIntensitiesSimObs];
    primaryRatiosSim = [primaryRatiosSim;primaryRatiosSimObs];
    
    % Recover observer parameters
    [calcParams,~] = findObserverParameters(testSpds,primarySpds,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize);
    recoveredParams = [recoveredParams;calcParams];
    
    % Calculate predicted test intensity and primary ratio based on
    % recovered parameters
    [~,~,testIntensitiesPredObs,primaryRatiosPredObs] = ...
        getMatchSeries(subjID,[calcParams 0],p1,p2,test,...
        'predicted','fieldSize',p.Results.fieldSize,...
        'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'saveResults',false);
    testIntensitiesPred = [testIntensitiesPred;testIntensitiesPredObs];
    primaryRatiosPred = [primaryRatiosPred;primaryRatiosPredObs];
    
    % Calculate root mean square error of the spectral sensitivities for
    % two sets of parameters, and for the sampled parameters compared to
    % the standard observer. 
    [coneErr(i),coneErr540(i)] = findConeSensitivityError(observerParams,...
        [calcParams 0],'age',p.Results.age,'fieldSize',p.Results.fieldSize);
    [coneStandardErr(i),coneStandardErr540(i)] = ...
        findConeSensitivityError(observerParams,zeros(1,9),'age',...
        p.Results.age,'fieldSize',p.Results.fieldSize);
end

% Output structs
error = struct();           % Cone spectra sensitivity error
error.coneErr = coneErr; 
error.coneStandardErr = coneStandardErr;
error.coneErr540 = coneErr540;
error.coneStandardErr540 = coneStandardErr540;

rayleighSettings = struct(); % Primary ratio and test intensity settings
rayleighSettings.testIntensitiesSim = testIntensitiesSim;
rayleighSettings.testIntensitiesPred = testIntensitiesPred;
rayleighSettings.primaryRatiosSim = primaryRatiosSim;
rayleighSettings.primaryRatiosPRed = primaryRatiosPred;

%% Plots
% Actual vs predicted plots for each parameter
coneParamNames = {'Lens Density','Macular Pigment Density',...
    'L Photopigment Density','M Photopigment Density',...
    'S photopigment density','L Lambda Max','M Lambda Max', 'S Lambda Max'};
for k = 1:length(paramsToVary)
    xVals = sampledParams(:,k);   % Predicted parameters
    yVals = recoveredParams(:,k); % Recovered params
    figure();
    hold on;
    if (k == 6) || (k==7) || (k==8)  % Lambda max shifts, in nm
        limits = [-4 4];
    else                             % Density shifts, in %
        limits = [-30 30]; 
    end
    xlim(limits);
    ylim(limits);
    l1 = plot(xVals,yVals,'b* ','MarkerSize',7,'LineWidth',1);
    l2 = refline(1,0);
    theTitle = sprintf('%s Predicted vs Actual',cell2mat(coneParamNames(k)));
    title(theTitle);
    xlabel('Sampled Parameters');
    ylabel('Recovered Parameters');
    legend([l1 l2],'Parameters','y=x');
end

% Comparative cone error plot - all wavelengths 
combinedConeErr = [coneErr coneStandardErr];
figure();
bar(combinedConeErr);
title('Cone Spectral Sensitivity Error');
legend('Sampled vs Recovered Params', 'Sampled vs Standard Params');
ylabel('Error'); 
xlabel('Observer'); 
if (length(method) == 9) && all(method=='predicted')
    ylim([0 0.005]);
else
    ylim([0 0.01]);
end 

% Comparative cone error plot - wavelengths over 540nm 
combinedConeErr540 = [coneErr540 coneStandardErr540];
figure();
bar(combinedConeErr540);
title('Cone Spectral Sensitivity Error Above 540nm');
legend('Sampled vs Recovered Params', 'Sampled vs Standard Params');
ylabel('Error'); 
xlabel('Observer');
ylim([0 0.01]);
if (length(method)==9) && all(method=='predicted')
    ylim([0 0.005]);
else
    ylim([0 0.01]);
end 

% Cone spectral sensitivities and generalized Pitt diagrams 
%
% For which observers are we making the plots?
if p.Results.makeAllPlots
    plottingInds = 1:nObservers;
else
    [~,bestObs] = min(coneErr);
    [~,worstObs] = max(coneErr);
    plottingInds = [bestObs worstObs];
end

for k = plottingInds
    % Cone spectral sensitivity plot
    S = [380 2 201];
    sampledObserver = genRayleighObserver('coneVec',...
        sampledParams(k,:),'S',S);
    recoveredObserver = genRayleighObserver('coneVec',...
        [recoveredParams(k,:) 0],'S',S);
    figure();
    hold on;
    l1 = plot(SToWls(S),sampledObserver.T_cones(1:2,:),'b-',...
        'LineWidth',2.5);
    l2 = plot(SToWls(S),recoveredObserver.T_cones(1:2,:),'r-',...
        'LineWidth',1.25);
    legend([l1(1) l2(1)],'sampled observer','recovered observer');
    theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g', k);
    title(theTitle);
    
    % Generalized Pitt diagram
    figure();
    hold on;
    l1 = plot(primaryRatiosSim(k,:),testIntensitiesSim(k,:),'b-o',...
        'LineWidth',2.5);
    l2 = plot(primaryRatiosPred(k,:),testIntensitiesPred(k,:),'r-o',...
        'LineWidth',1.25);
    theTitle = sprintf('Generalized Pitt Diagram, observer %g',k);
    title(theTitle);
    xlabel('Primary Ratio');
    ylabel('Test Intensity');
    legend([l1 l2],'Simulated','Predicted');
    
    % Add wavelength labels to Pitt diagram (assume test light is shifted)
    labels = test';
    labelsText = cellstr(num2str(labels));
    dx =  -0.01;
    dy = 0.03;
    text(primaryRatiosSim(k,:)+dx,testIntensitiesSim(k,:)+dy,...
        labelsText);
end

% If individual plots were made for best/worst observers, modify plot
% titles and highlight points
if ~p.Results.makeAllPlots
    figure(length(paramsToVary)+3);
    theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g (best)',bestObs);
    title(theTitle); 
    
    figure(length(paramsToVary)+4);
    theTitle = sprintf('Generalized Pitt Diagram, observer %g (best)',bestObs);   
    title(theTitle); 
    
    figure(length(paramsToVary)+5);    
    theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g (worst)',worstObs);
    title(theTitle); 

    figure(length(paramsToVary)+6);
    theTitle = sprintf('Generalized Pitt Diagram, observer %g (worst)',worstObs);   
    title(theTitle); 
    
    for i = 1:length(paramsToVary)
        figure(i) 
        plot(sampledParams(bestObs,i),recoveredParams(bestObs,i),'gs',...
            'MarkerSize',10,'LineWidth',2);
        plot(sampledParams(worstObs,i),recoveredParams(worstObs,i),'rs',...
            'MarkerSize',8,'LineWidth',1.5);
        legend('parameter','y = x','best observer','worst observer');
    end 
end
end 