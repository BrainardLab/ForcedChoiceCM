function [coneAvgErr,matchAvgErr] = sampleRayleighMatch(subjID,...
    nObservers,baseParams,paramsToVary,p1,p2,test,method,varargin)
% Tests the OLRayleighMatch simulation and makes associated sampling plots.

% Syntax:
%   sampledRayleighMatch(subjID,nObservers,baseParams,paramsToVary,
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
%                         outlined above and saves as pdfs. When false,
%                         makes no plots. Default is true.
%    'makeAllObserverPlots' -Logical. When true, plots the cone spectral
%                            sensitivities and generalized Pitt diagrams
%                            for each observer, not just the best and
%                            worst. Default is false.
%    'fieldSize          -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'adjustmentLength'  -Integer defining the size of the lights array
%                         available for OLRayleighMatch. Default is 3201.
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
%    'dmac0'             -Logical. If true, the parameter search constrains
%                         macular pigment density to 0. Default is false.
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
%   07/23/20  dce       Changed plotting
%   07/24/20  dce       Added option to lock macular pigment 

% Close stray figures
close all;

% Parse input
p = inputParser;
p.addParameter('plotResults',true,@(x)(islogical(x)));
p.addParameter('makeAllObserverPlots',false,@(x)(islogical(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',3201,@(x)(isnumeric(x)));
p.addParameter('monochromatic',false,@(x)(islogical(x)));
p.addParameter('nominal',false,@(x)(islogical(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',false,@(x)(islogical(x)));
p.addParameter('dmac0',false,@(x)(islogical(x)));
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
optErrs = 0;                           % Counter for optimization errors

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
        'saveResults',false,'nominal',p.Results.nominal,'adjustmentLength',...
        p.Results.adjustmentLength);
    testIntensitiesSim = [testIntensitiesSim;testIntensitiesSimObs];
    primaryRatiosSim = [primaryRatiosSim;primaryRatiosSimObs];
    
    % Record whether getMatchSeries made successful matches in all cases
    nAttemptedMatches = length(test);    % Number of attempted matches
    [~,nMatchSpds] = size(testSpds);     % Number of successful matches
    
    % Recover observer parameters and associated match error
    try 
        [calcParams,calcErr] = findObserverParameters(testSpds,primarySpds,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'restrictBySd',p.Results.restrictBySd,'dlens0',p.Results.dlens0,...
        'LMEqualOD',p.Results.LMEqualOD,'dmac0',p.Results.dmac0,....
        'initialParams',baseParams,'errScalar',matchErrScalar);
    catch 
        warning('fmincon searched for impossible cone parameters. Currently redoing search');
        optErrs = optErrs+1;
        [calcParams,calcErr] = findObserverParameters(testSpds,primarySpds,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'restrictBySd',p.Results.restrictBySd,'dlens0',p.Results.dlens0,...
        'LMEqualOD',p.Results.LMEqualOD,'dmac0',p.Results.dmac0,....
        'initialParams',baseParams,'errScalar',matchErrScalar);
    end 
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
save(fullfile(outputDir,'paramsSearchData.mat'));

%% Plots - produces two pdfs with subplots
if p.Results.plotResults
    %% Figure 1 - error plots
    % How many subplots are we making?
    if p.Results.makeAllObserverPlots
        plottingInds = 1:nObservers;
    else  % Find the indices of the best and worst observers
        [~,bestObs] = min(coneErr);
        [~,worstObs] = max(coneErr);
        plottingInds = [bestObs worstObs];
    end
    nErrPlots = 2*length(plottingInds)+2; % Number of subplots
    nCols1 = 2;                           % Number of columns
    nRows1 = ceil(nErrPlots/nCols1);      % Number of rows
    
    % Set up the base figure
    errFig = figure(1);
    set(errFig,'Color',[1 1 1],'Position',[10 10 1400 800]);
    hold on;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', nRows1, ...
        'colsNum', nCols1, ...
        'heightMargin',  0.1, ...
        'widthMargin',    0.1, ...
        'leftMargin',     0.07, ...
        'rightMargin',    0.04, ...
        'bottomMargin',   0.1, ...
        'topMargin',      0.04);
    
    % Subplot 1 - cone spectral sensitivity error
    subplot('Position', subplotPosVectors(1,1).v);
    hold on;
    bar([coneErr coneStandardErr]);
    title('Cone Spectral Sensitivity Error');
    legend('Sampled vs Recovered Params', 'Sampled vs Base Params');
    ylabel('Error');
    xlabel('Observer');
    ylim([0 0.02]);
    
    % Subplot 2 - match error
    subplot('Position', subplotPosVectors(1,2).v);
    hold on;
    bar([matchSampledErr matchErr matchStandardErr]);
    title('Match Spectral Sensitivity Error');
    legend('Sampled Params', 'Recovered Params', 'Base Params');
    ylabel('Error');
    xlabel('Observer');
    ylim([0 0.2]);
    
    % Additional subplots - cone spectra and generalized Pitt diagrams
    S = [380 2 201];
    for k = 1:length(plottingInds)
        row = k+1;     % Row of current subplot
        % Cone plot
        subplot('Position', subplotPosVectors(row,1).v);
        hold on;
        % Find and plot observer cone fundamentals
        sampledObserver = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'coneVec',...
            sampledParams(plottingInds(k),:));
        recoveredObserver = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'coneVec',...
            [recoveredParams(plottingInds(k),:) 0]);
        l1 = plot(SToWls(S),sampledObserver.T_cones(1:2,:),'b-',...
            'LineWidth',2.5);
        l2 = plot(SToWls(S),recoveredObserver.T_cones(1:2,:),'r-',...
            'LineWidth',1.25);
        % Clean up plot
        legend([l1(1) l2(1)],'sampled observer','recovered observer');
        theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g',...
            plottingInds(k));
        title(theTitle);
        xlabel('Wavelength (nm)');
        ylabel('Power');
        
        % Generalized Pitt diagram - plot simulated and predicted match
        % data
        subplot('Position', subplotPosVectors(row,2).v);
        hold on;
        l1 = plot(primaryRatiosSim(plottingInds(k),:),...
            testIntensitiesSim(plottingInds(k),:),'b-o','LineWidth',2.5);
        l2 = plot(primaryRatiosPred(plottingInds(k),:),...
            testIntensitiesPred(plottingInds(k),:),'r-o','LineWidth',1.25);
        % Clean up plot
        theTitle = sprintf('Generalized Pitt Diagram, observer %g',plottingInds(k));
        title(theTitle);
        xlabel('Primary Ratio');
        ylabel('Test Intensity');
        lgd = legend([l1 l2],'Simulated','Predicted');
        lgd.Location = 'northwest';
        % Add text labels for wavelength
        labels = cellstr(num2str(test'));
        dx =  -0.01;   % x offset
        dy = 0.005;     % y offset
%         text(primaryRatiosSim(plottingInds(k),:)+dx...
%             ,testIntensitiesSim(plottingInds(k),:)+dy,labels);
    end
    % Edit titles if plots were only made for best and worst observers
    if ~p.Results.makeAllObserverPlots
        subplot('Position', subplotPosVectors(2,1).v);
        theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g (best)',bestObs);
        title(theTitle);
        
        subplot('Position', subplotPosVectors(2,2).v);
        theTitle = sprintf('Generalized Pitt Diagram, observer %g (best)',bestObs);
        title(theTitle);
        
        subplot('Position', subplotPosVectors(3,1).v);
        theTitle = sprintf('L and M Cone Spectral Sensitivities, observer %g (worst)',worstObs);
        title(theTitle);
        
        subplot('Position', subplotPosVectors(3,2).v);
        theTitle = sprintf('Generalized Pitt Diagram, observer %g (worst)',worstObs);
        title(theTitle);
    end
    % Save plot
    sgtitle([subjID ' Error'],'Interpreter', 'none');
    NicePlot.exportFigToPDF(fullfile(outputDir,'errPlots'),...
        errFig,300);
    %% Figure 2 - parameter recovery
    paramsFig = figure(2);
    set(paramsFig,'Color',[1 1 1],'Position',[10 10 1700 800]);
    hold on;
    nCols2 = 4;
    nRows2 = ceil(length(paramsToVary)/nCols2);
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', nRows2, ...
        'colsNum', nCols2, ...
        'heightMargin',  0.07, ...
        'widthMargin',    0.07, ...
        'leftMargin',     0.04, ...
        'rightMargin',    0.04, ...
        'bottomMargin',   0.07, ...
        'topMargin',      0.04);
    coneParamNames = {'Lens Density','Macular Pigment Density',...
        'L Photopigment Density','M Photopigment Density',...
        'S photopigment density','L Lambda Max','M Lambda Max', 'S Lambda Max'};
    for ii = 1:length(paramsToVary)
        % Make a subplot in the correct position
        row = ceil(ii/nCols2);
        col = mod(ii,nCols2);
        if col == 0
            col = nCols2;
        end
        subplot('Position', subplotPosVectors(row,col).v);
        hold on;
        
        % Define axis limits
        if (ii==6) || (ii==7) || (ii==8)  % Lambda max shifts, in nm
            limits = [-5 5];
        else                             % Density shifts, in percent
            limits = [-40 40];
        end
        xlim(limits);
        ylim(limits);
        axis('square');
        
        % Plot data
        xVals = sampledParams(:,ii);   % Predicted parameters
        yVals = recoveredParams(:,ii); % Recovered params
        l1 = plot(xVals,yVals,'b* ','MarkerSize',7,'LineWidth',1);
        l2 = refline(1,0);
        
        % Titles and labels
        theTitle = sprintf('%s Predicted vs Actual',cell2mat(coneParamNames(ii)));
        title(theTitle);
        xlabel('Sampled Parameters');
        ylabel('Recovered Parameters');
        lgd = legend([l1 l2],'Parameters','y=x');
        lgd.Location = 'northwest';
        
        % If error plots were made for the best and worst observers, highlight
        % these observers on the parameter plots
        if ~p.Results.makeAllObserverPlots
            plot(sampledParams(bestObs,ii),recoveredParams(bestObs,ii),'gs',...
                'MarkerSize',10,'LineWidth',2);
            plot(sampledParams(worstObs,ii),recoveredParams(worstObs,ii),'rs',...
                'MarkerSize',8,'LineWidth',1.5);
            legend('parameter','y = x','best observer','worst observer');
        end
        
        % Save figure
        sgtitle([subjID ' Parameters'],'Interpreter', 'none');
        NicePlot.exportFigToPDF(fullfile(outputDir,'paramPlots'),...
            paramsFig,300);
    end
end
end