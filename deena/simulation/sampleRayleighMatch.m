function [coneAvgErr,matchAvgErr] = sampleRayleighMatch(subjID,...
    nObservers,baseParams,paramsToVary,p1,p2,test,method,varargin)
% Tests the OLRayleighMatch simulation and makes associated sampling plots.

% Syntax:
%   plotSampledRayleighMatch(subjID,nObservers,baseParams,paramsToVary,
%   p1,p2,test,method)
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
%    The program produces five types of plots. The first is a plot of
%    sampled vs predicted values for each parameter, aggregated across all
%    observers. The second is a bar plot of the difference between the cone
%    sensitivities of the sampled and recovered parameters for each
%    observer (root mean square error), compared to the difference for
%    the sampled parameters and the standard observer. The third is a bar
%    plot of the match error (Euclidian distance of the opponent contrast,
%    taken as rms error across all matches) for the sampled, recovered, and
%    standard observers. The fourth is a comparative plot of cone spectral
%    sensitivities based on the sampled and recovered parameters for each
%    observer. The fifth, a "generalized Pitt diagram," plots primary ratio
%    and test intensity settings that were found in the simulation for each
%    observer, comparing them to the settings predicted analyticall from
%    the recovered parameters.
%
%    By default, the fourth annd fifth plot types are produced for the best
%    and the worst observers (identified based on how different their cone
%    spectral sensitivities are). However, there is also an option to
%    produce the plots for every observer.
%
% Inputs:
%    subjID         -Character vector of subject ID
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
%    coneAvgErr     -Average difference between sampled and recovered cone
%                    spectral sensitivities, across all observers sampled.
%    matchAvgErr    -Average error associated with the simulated
%                    matches when using the recovered cone parameters.
%
% Optional key-value pairs:
%    'plotResults'       -Logical. When true, makes the five plot types
%                         outlined above. When false, makes no plots.
%                         Default is true.
%    'makeAllPlots'      -Logical. When true, plots the cone spectral
%                         sensitivities and generalized Pitt diagrams for
%                         each observer, not just the best and worst.
%                         Default is false.
%    'fieldSize          -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'monochromatic'     -When using analytic matches, logical indicating
%                         to use monochromatic spds insted of those
%                         generated by the OneLight. Default is false.
%    'nominal'           -Logical indicating to run the simulation with
%                         nominal, rather than predicted, spds. Default
%                         is false.
%    'nReversals'        -Number of reversals the simulated observer must
%                         make before changing step size. Enter as a 2-
%                         element vector - the first element is the number
%                         of reversals for intermediate step sizes, the
%                         second is the number needed for the smallest step
%                         size. Default is [1 4].
%    'nObserverMatches'  -Number of matches to make with each simulated set
%                         of lights. Default is 1.
%    'nBelowThreshold'   -When using a simulated observer with
%                         threshold matching, number of pairs below
%                         threshold required before recording a match.
%                         Default is 1.
%    'thresholdScaleFactor' -When using a simulated observer with
%                            threshold matching, scale factor for matching
%                            threshold. Default is 0.5.
%    'LMEqualOD'         -Logical. If true, the parameter search constrains
%                         L and M cone optical densities to be equal.
%                         Default is false.
%    'dlens0'            -Logical. If true, the parameter search constrains
%                         the lens pigment density to 0. Default is false.
%    'restrictBySd'      -Logical. If true, the parameter search restricts
%                         all params to within three standard deviations of
%                         their means. Default is true.

% History:
%   07/06/20  dce       Wrote it.
%   07/10/20  dce       Edited plotting and added rms error output.
%   07/16/20  dce       Got rid of 540nm cone error plots
%   07/17/20  dce       Added saving
%   07/21/20  dce       Style edits, added param recovery key-value pairs
%   07/22/20  dce       Added error outputs

% Close stray figures
close all;

% Parse input
p = inputParser;
p.addParameter('plotResults',true,@(x)(islogical(x)));
p.addParameter('makeAllPlots',false,@(x)(islogical(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('monochromatic',false,@(x)(islogical(x)));
p.addParameter('nominal',false,@(x)(islogical(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',false,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Standard deviation of parameters (Asano 2015)
sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];
matchErrScalar = 100;    % Scale factor to improve match error search 

% Base observer, used for comparison
stdObs = genRayleighObserver('age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'coneVec',baseParams);

% Create directory for saving results
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
else
    error('Specified directory already exists')
end

% Data-storing arrays. Each row holds a different observer's params.
sampledParams = [];           % Sampled cone parameters
recoveredParams = [];         % Cone parameters recovered by simulation
testIntensitiesSim = [];      % Calculated test intensities
testIntensitiesPred = [];     % Test intensities predicted from recovered params
primaryRatiosSim = [];        % Calculated primary ratios
primaryRatiosPred = [];       % Primary ratios predicted from recovered params
coneErr = zeros(nObservers,1);         % Error between sampled and recovered cones
coneStandardErr = zeros(nObservers,1); % Error between sampled and standard cones
matchErr = zeros(nObservers,1);        % Match error when using recovered parameters
matchSampledErr = zeros(nObservers,1); % Match error when using sampled parameters
matchStandardErr = zeros(nObservers,1);% Match error when using base parameters

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
    
    % Make a series of Rayleigh matches for the observer
    testingID = 'test_series';
    [testSpds,primarySpds,testIntensitiesSimObs,primaryRatiosSimObs] = ...
        getMatchSeries(testingID,observerParams,p1,p2,test,method,'fieldSize',...
        p.Results.fieldSize,'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'nObserverMatches',...
        p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,'thresholdScaleFactor',...
        p.Results.thresholdScaleFactor,'rayleighPlots',false,...
        'saveResults',false,'nominal',p.Results.nominal);
    testIntensitiesSim = [testIntensitiesSim;testIntensitiesSimObs];
    primaryRatiosSim = [primaryRatiosSim;primaryRatiosSimObs];
    
    % Recover observer parameters and associated match error
    [calcParams,calcErr] = findObserverParameters(testSpds,primarySpds,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'restrictBySd',p.Results.restrictBySd,'dlens0',p.Results.dlens0,...
        'LMEqualOD',p.Results.LMEqualOD,'initialParams',baseParams,...
        'errScalar',matchErrScalar);
    recoveredParams = [recoveredParams;calcParams];
    matchErr(i) = calcErr;
    
    % Recover match error for the standard observer and the sampled
    % observer
    matchStandardErr(i) = findMatchError(zeros(1,8),stdObs,testSpds,...
        primarySpds,'errScalar',matchErrScalar)/matchErrScalar;
    matchSampledErr(i) = findMatchError(observerParams(1:8),...
        genRayleighObserver('coneVec',observerParams,'age',p.Results.age,...
        'fieldSize',p.Results.fieldSize),testSpds,primarySpds,...
        'errScalar',matchErrScalar)/matchErrScalar;
    
    % Calculate predicted test intensity and primary ratio based on
    % recovered parameters. This will be compared to the actual values
    % found by the simulation.
    [~,~,testIntensitiesPredObs,primaryRatiosPredObs] = ...
        getMatchSeries(testingID,[calcParams 0],p1,p2,test,...
        'predicted','fieldSize',p.Results.fieldSize,...
        'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'saveResults',false,...
        'nominal',p.Results.nominal);
    testIntensitiesPred = [testIntensitiesPred;testIntensitiesPredObs];
    primaryRatiosPred = [primaryRatiosPred;primaryRatiosPredObs];
    
    % Calculate root mean square error of the spectral sensitivities for
    % the two sets of parameters, and for the sampled parameters compared
    % to the base observer.
    [coneErr(i)] = findConeSensitivityError(observerParams,...
        [calcParams 0],'age',p.Results.age,'fieldSize',p.Results.fieldSize);
    [coneStandardErr(i)] = findConeSensitivityError(observerParams,...
        baseParams,'age',p.Results.age,'fieldSize',p.Results.fieldSize);
    save(fullfile(outputDir,'paramsSearchData.mat'));
end

% Compute average error measures 
coneAvgErr = mean(coneErr);   % Average cone spectral sensitivity error
matchAvgErr = mean(matchErr); % Average match error with recovered params

%% Plots
if p.Results.plotResults
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
            limits = [-5 5];
        else                             % Density shifts, in percent
            limits = [-40 40];
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
    
    % Comparative cone error plot
    figure();
    bar([coneErr coneStandardErr]);
    title('Cone Spectral Sensitivity Error');
    legend('Sampled vs Recovered Params', 'Sampled vs Base Params');
    ylabel('Error');
    xlabel('Observer');
    ylim([0 0.02]);
    
    % Comparative match error plot
    figure();
    bar([matchSampledErr matchErr matchStandardErr]);
    title('Match Spectral Sensitivity Error');
    legend('Sampled Params', 'Recovered Params', 'Base Params');
    ylabel('Error');
    xlabel('Observer');
    ylim([0 0.2]);
    
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
        sampledObserver = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'coneVec',...
            sampledParams(k,:));
        recoveredObserver = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'coneVec',...
            [recoveredParams(k,:) 0]);
        figure();
        hold on;
        l1 = plot(SToWls(S),sampledObserver.T_cones(1:2,:),'b-',...
            'LineWidth',2.5);
        l2 = plot(SToWls(S),recoveredObserver.T_cones(1:2,:),'r-',...
            'LineWidth',1.25);
        legend([l1(1) l2(1)],'sampled observer','recovered observer');
        theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g', k);
        title(theTitle);
        xlabel('Wavelength (nm)');
        ylabel('Power');
        
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
        labels = cellstr(num2str(test'));
        dx =  -0.01;
        dy = 0.03;
        text(primaryRatiosSim(k,:)+dx,testIntensitiesSim(k,:)+dy,labels);
    end
    
    % If individual plots were made for best/worst observers, modify plot
    % titles and highlight points
    if ~p.Results.makeAllPlots
        for i = 1:length(paramsToVary)
            figure(i)
            plot(sampledParams(bestObs,i),recoveredParams(bestObs,i),'gs',...
                'MarkerSize',10,'LineWidth',2);
            plot(sampledParams(worstObs,i),recoveredParams(worstObs,i),'rs',...
                'MarkerSize',8,'LineWidth',1.5);
            legend('parameter','y = x','best observer','worst observer');
        end
        
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
    end
    
    % Save figures
    figs =  findobj('type','figure');
    for i = 1:length(figs)
        figure(i)
        print(fullfile(outputDir,[subjID '_' num2str(i)]),'-dtiff');
    end
end
end