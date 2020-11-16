function [coneAvgErr,matchAvgErr,coneAvgStdErr,matchAvgStdErr,...
    matchAvgSampledErr,sampledParams,recoveredParams]...
    = sampleRayleighMatch(subjID,nObservers,baseParams,paramsToVary,...
    opponentParams,p1,p2,test,method,varargin)
% Tests the OLRayleighMatch simulation and makes associated sampling plots.
%
% Syntax:
%   sampleRayleighMatch(subjID,nObservers,baseParams,paramsToVary,
%   opponentParams,p1,p2,test,method)
%
% Description:
%    Takes in a vector of cone individual difference parameters as a
%    starting point, as well as information about which parameters to vary.
%    Creates a series of simulated observers by sampling the chosen
%    parameters based on their standard deviations. For each observer, runs
%    a Rayleigh matching simulation with a set of primary/test lights
%    and uses the match results to recover the individual difference
%    parameters (findObserverParameters). Then, conducts a few basic
%    analyses to validate how similar the recovered parameters are to the
%    sampled (true) parameters.
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
%    baseParams     -Eight-element numeric vector of individual difference
%                    parameters (see ObserverVecToParams for description).
%                    The values are used as means for observer sampling.
%    paramsToVary   -Eight-element numeric vector of ones and zeros
%                    indicating which individual difference parameters
%                    should be varied (the noise parameter is excluded).
%                    Parameters set to 1 will be sampled around their
%                    standard deviation, while parameters set to 0 will
%                    stay at the values specified in baseParams.
%    opponentParams -Four-element vector of opponent contrast weightings.
%                    (see genRayleighObserver).
%    p1             -Integer or numeric vector of desired wavelengths
%                    for the first primary.
%    p2             -Integer or numeric vector of desired wavelengths
%                    for the second primary.
%    test           -Integer or numeric vector of desired wavelengths
%                    for the test light.
%    method         -Character vector for match method. Choose either
%                    'predicted', 'forcedChoice','threshold', or
%                    'bestAvailable'.
%
% Outputs:
%    coneAvgErr     -Average difference between sampled and recovered cone
%                    spectral sensitivities, across all observers sampled.
%    matchAvgErr    -Average error associated with the simulated
%                    matches when using the recovered cone parameters.
%    coneAvgStdErr  -Average difference between sampled and standard cone
%                    spectral sensitivities, across all observers sampled.
%    matchAvgStdErr -Average error associated with the simulated
%                    matches when using the standard cone parameters.
%    matchAvgSampledErr -Average error associated with the simulated
%                        matches when using the sampled cone parameters.
%    sampledParams  -nObservers x length(baseParams) matrix of sampled
%                    observer parameters.
%    recoveredParams-nObservers x length(baseParams) matrix of recovered
%                    observer parameters.
%
% Optional key-value pairs:
%    'makeAllObserverPlots' -Logical. When true, plots the cone spectral
%                            sensitivities and generalized Pitt diagrams
%                            for each observer, not just the best and
%                            worst. Default is false.
%    'makeNoObserverPlots'  -Logical. When true, no plots are made. Default
%                            is false.
%    'fieldSize'         -Integer field size, in degrees. Default is 2.
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
%    'monochromatic'     -Logical indicating to use monochromatic spds
%                         instead of those generated by the OneLight.
%                         Default is false.
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
%    'noiseScaleFactor'  -Number >=0 which determines which scalar
%                         multiple of the opponent noise SD should be
%                         used as the observer noise SD. Default is 0.
%    'LMEqualOD'         -Logical. If true, the parameter search constrains
%                         L and M cone optical densities to be equal.
%                         Default is false.
%    'dlens0'            -Logical. If true, the parameter search constrains
%                         the lens pigment density to 0. Default is true.
%    'dmac0'             -Logical. If true, the parameter search constrains
%                         macular pigment density to 0. Default is true.
%    'restrictBySd'      -Logical. If true, the parameter search restricts
%                         all params to within two standard deviations of
%                         their means. Default is true.
%    'sampledObservers'  -nObservers x length(baseParams)array of
%                         previously-sampled observer parameters (useful
%                         for when observers are used across multiple
%                         experimental conditions). Default is [].
%    'errWls'            -numeric vector of test wavelengths. If it is
%                         nonempty, average match errors are calculated
%                         only for these wavelengths. Default is [].
%    'plotTitle'         -Character vector for titling plots. Default is
%                         [], in which case a default file naming system is
%                         used.
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on 
%                         stimulus parameters. Each row represents a given  
%                         test wavelength and the limits which are associated 
%                         with it. The columns are arranged as follows: 
%                         [test wl, min lambda, max lambda, min test 
%                         intensity, max test intensity]. Default is [].
%    'lambdaRef'         - Number between 0 and 1 indicating which value 
%                          of lambda to use for a reference primary when
%                          calculating simulated opponent contrasts. Must 
%                          be a member of p1Scales. When empty, opponent 
%                          contrasts are not computed relative to a 
%                          reference. Default is []. 

% History:
%   07/06/20  dce       Wrote it.
%   07/10/20  dce       Edited plotting and added rms error output.
%   07/16/20  dce       Got rid of 540nm cone error plots
%   07/17/20  dce       Added saving
%   07/21/20  dce       Style edits, added param recovery key-value pairs
%   07/22/20  dce       Added error outputs
%   07/23/20  dce       Changed plotting
%   07/24/20  dce       Added option to lock macular pigment
%   07/29/20  dce       Moved parameter sampling to a separate function
%   07/31/20  dce       Changed plot titling and filename conventions
%   08/03/20  dce       Added more outputs (standard and sampled)
%   08/09/20  dce       Changed generalized Pitt diagram to include best
%                       available matches, not analytic
%   10/28/20  dce       Added sampled and recovered parameters as outputs
%   11/01/20  dce       Added wavelength sampling as key-value pair
%   11/15/20  dce       Added option to restrict stimuli

% Close stray figures
close all;

% Parse input
p = inputParser;
p.addParameter('makeAllObserverPlots',false,@(x)(islogical(x)));
p.addParameter('makeNoObserverPlots',false,@(x)(islogical(x)));
p.addParameter('plotTitle',[],@(x)(ischar(x)));
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
p.addParameter('noiseScaleFactor',0,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('errWls',[],@(x)(isnumeric(x)));
p.addParameter('sampledObservers',[],@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.addParameter('lambdaRef',[],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Base observer, used for comparison
stdObs = genRayleighObserver('age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'coneVec',baseParams,'opponentParams',...
    opponentParams,'S',p.Results.S);

% Create directory for saving results
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
else
    error('Specified directory already exists')
end

% Get sampled parameters. Each row represents a different observer's params
if isempty(p.Results.sampledObservers)
    sampledParams = sampleRayleighObservers(nObservers,baseParams,...
        paramsToVary);
else
    sampledParams = p.Results.sampledObservers;
    [r,c] = size(sampledParams);
    if r~=nObservers || c~=length(baseParams)
        error('Provided observer params have incorrect dimensions');
    end
end

% Initial parameters related to error calculation
matchErrScalar = 100;    % Scale factor to improve match error search
if ~isempty(p.Results.errWls)  % Calculate average error for specific wls
    selectionArr = zeros(1,length(test));
    for i = 1:length(p.Results.errWls)
        selectionArr = selectionArr + (test==p.Results.errWls(i));
    end
else                           % Calculate average error using all wls
    selectionArr = ones(1,length(test));
end
selectionArr = logical(selectionArr);

% Data-storing arrays. Each row holds a different observer's params.
recoveredParams = [];         % Cone parameters recovered by simulation
testIntensitiesSim = [];      % Calculated test intensities
testIntensitiesRecPred = [];     % Test intensities predicted from recovered params
testIntensitiesSimPred = [];     % Test intensities predicted from simulated params
primaryRatiosSim = [];            % Calculated primary ratios
primaryRatiosRecPred = [];       % Primary ratios predicted from recovered params
primaryRatiosSimPred = [];       % Primary ratios predicted from simulated params
coneErr = zeros(nObservers,1);         % Error between sampled and recovered cones
coneStandardErr = zeros(nObservers,1); % Error between sampled and standard cones
matchErr = zeros(nObservers,1);        % Match error when using recovered parameters
matchSampledErr = zeros(nObservers,1); % Match error when using sampled parameters
matchStandardErr = zeros(nObservers,1);% Match error when using base parameters

% For each observer: sample parameters, make matches, and use matches to
for i = 1:nObservers
    % Make a series of Rayleigh matches for the observer
    testingID = [subjID '_' num2str(i)];
    [testSpds,primarySpds,testIntensitiesSimObs,primaryRatiosSimObs] = ...
        getMatchSeries(testingID,sampledParams(i,:),opponentParams,p1,p2,...
        test,method,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
        'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,...
        'testScale',p.Results.testScale,'monochromatic',...
        p.Results.monochromatic,'nObserverMatches',...
        p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,'thresholdScaleFactor',...
        p.Results.thresholdScaleFactor,'rayleighPlots',false,...
        'saveResults',false,'nominal',p.Results.nominal,'adjustmentLength',...
        p.Results.adjustmentLength,'noiseScaleFactor',...
        p.Results.noiseScaleFactor,'averageSpds',true,'sPredicted',...
        p.Results.S,'stimLimits',p.Results.stimLimits,'lambdaRef',...
        p.Results.lambdaRef);
    testIntensitiesSim = [testIntensitiesSim;testIntensitiesSimObs];
    primaryRatiosSim = [primaryRatiosSim;primaryRatiosSimObs];
    testSpdsSelected = testSpds(:,selectionArr);
    primarySpdsSelected = primarySpds(:,selectionArr);
    
    % Recover observer parameters and associated match error
    [calcParams] = findObserverParameters(testSpds,primarySpds,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'restrictBySd',p.Results.restrictBySd,'dlens0',p.Results.dlens0,...
        'LMEqualOD',p.Results.LMEqualOD,'dmac0',p.Results.dmac0,....
        'initialConeParams',baseParams,'errScalar',matchErrScalar,...
        'opponentParams',opponentParams,'S',p.Results.S);
    recoveredParams = [recoveredParams;calcParams];
    
    % Find match error for the recovered, standard, and sampled observers
    matchErr(i) = findMatchError(calcParams,stdObs,testSpdsSelected,...
        primarySpdsSelected,'errScalar',matchErrScalar,'S',p.Results.S...
        )/matchErrScalar;
    matchStandardErr(i) = findMatchError(zeros(1,8),stdObs,testSpdsSelected,...
        primarySpdsSelected,'errScalar',matchErrScalar,...
        'S',p.Results.S)/matchErrScalar;
    matchSampledErr(i) = findMatchError(sampledParams(i,:),stdObs,...
        testSpdsSelected,primarySpdsSelected,'errScalar',matchErrScalar,...
        'S',p.Results.S)/matchErrScalar;
    
    % Calculate predicted test intensity and primary ratio based on
    % recovered parameters. This will be compared to the actual values
    % found by the simulation.
    [~,~,testIntensitiesPredObs,primaryRatiosPredObs] = ...
        getMatchSeries(testingID,calcParams,opponentParams,p1,p2,test,...
        'bestAvailable','fieldSize',p.Results.fieldSize,...
        'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'saveResults',false,...
        'nominal',p.Results.nominal,'sPredicted',p.Results.S,'lambdaRef',...
        p.Results.lambdaRef);
    testIntensitiesRecPred = [testIntensitiesRecPred;testIntensitiesPredObs];
    primaryRatiosRecPred = [primaryRatiosRecPred;primaryRatiosPredObs];
    
    % Calculate predicted test intensity and primary ratio based on
    % simulated parameters.
    [~,~,testIntensitiesPredObs2,primaryRatiosPredObs2] = ...
        getMatchSeries(testingID,sampledParams(i,:),opponentParams,p1,p2,test,...
        'bestAvailable','fieldSize',p.Results.fieldSize,...
        'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'monochromatic',p.Results.monochromatic,'saveResults',false,...
        'nominal',p.Results.nominal,'sPredicted',p.Results.S,'lambdaRef',...
        p.Results.lambdaRef);
    testIntensitiesSimPred = [testIntensitiesSimPred;testIntensitiesPredObs2];
    primaryRatiosSimPred = [primaryRatiosSimPred;primaryRatiosPredObs2];
    
    % Calculate root mean square error of the spectral sensitivities for
    % the two sets of parameters, and for the sampled parameters compared
    % to the base observer.
    [coneErr(i)] = findConeSensitivityError(sampledParams(i,:),...
        calcParams,'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'opponentParams',opponentParams,'S',p.Results.S);
    [coneStandardErr(i)] = findConeSensitivityError(sampledParams(i,:),...
        baseParams,'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'opponentParams',opponentParams,'S',p.Results.S);
    if i~=nObservers
        save(fullfile(outputDir,[subjID '_paramsSearchData.mat']));
    end
end

% Compute average error measures
coneAvgErr = mean(coneErr);
coneAvgStdErr = mean(coneStandardErr);
matchAvgErr = mean(matchErr);
matchAvgStdErr = mean(matchStandardErr);
matchAvgSampledErr = mean(matchSampledErr);
save(fullfile(outputDir,[subjID '_paramsSearchData.mat']));

%% Plots - produces two pdfs with subplots
%% Figure 1 - error plots
% How many subplots are we making?
if ~p.Results.makeNoObserverPlots
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
        'topMargin',      0.1);
    
    % Subplot 1 - cone spectral sensitivity error
    subplot('Position', subplotPosVectors(1,1).v);
    hold on;
    bar([coneErr coneStandardErr]);
    title('Cone Spectral Sensitivity Error');
    legend('Simulated vs Recovered Params', 'Simulated vs Standard Params');
    ylabel('Error');
    xlabel('Observer');
    ylim([0 0.02]);
    
    % Subplot 2 - match error
    subplot('Position', subplotPosVectors(1,2).v);
    hold on;
    bar([matchSampledErr matchErr matchStandardErr]);
    title('Match Spectral Sensitivity Error');
    legend('Simulated Parameters','Recovered Parameters','Standard Parameters');
    ylabel('Error');
    xlabel('Observer');
    ylim([0 0.2]);
    
    % Additional subplots - cone spectra and generalized Pitt diagrams
    for k = 1:length(plottingInds)
        row = k+1;     % Row of current subplot
        % Cone plot
        subplot('Position', subplotPosVectors(row,1).v);
        hold on;
        % Find and plot observer cone fundamentals
        sampledObserver = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'coneVec',...
            sampledParams(plottingInds(k),:),'opponentParams',opponentParams,...
            'S',p.Results.S);
        recoveredObserver = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'coneVec',...
            recoveredParams(plottingInds(k),:),'opponentParams',opponentParams,...
            'S',p.Results.S);
        l1 = plot(SToWls(p.Results.S),sampledObserver.T_cones(1:2,:),'b-',...
            'LineWidth',2.5);
        l2 = plot(SToWls(p.Results.S),recoveredObserver.T_cones(1:2,:),'r-',...
            'LineWidth',1.25);
        % Clean up plot
        legend([l1(1) l2(1)],'Simulated Observer','Recovered Observer');
        theTitle = sprintf('L and M Cone Spectral Sensitivities, Observer %g',...
            plottingInds(k));
        title(theTitle);
        xlabel('Wavelength (nm)');
        ylabel('Power');
        
        % Generalized Pitt diagram - plot simulated and predicted match
        % data
        subplot('Position', subplotPosVectors(row,2).v);
        hold on;
        l1 = plot(primaryRatiosRecPred(plottingInds(k),:),...
            testIntensitiesRecPred(plottingInds(k),:),'r o','MarkerSize',8,...
            'MarkerFaceColor','r');
        l2 = plot(primaryRatiosSimPred(plottingInds(k),:),...
            testIntensitiesSimPred(plottingInds(k),:),'g o','MarkerSize',6,...
            'MarkerFaceColor','g');
        l3 = plot(primaryRatiosSim(plottingInds(k),:),...
            testIntensitiesSim(plottingInds(k),:),'b o','MarkerSize',4,...
            'MarkerFaceColor','b');
        % Clean up plot
        theTitle = sprintf('Generalized Pitt Diagram, Observer %g',plottingInds(k));
        title(theTitle);
        xlabel('Primary Ratio');
        ylabel('Test Intensity');
        lgd = legend([l1 l2 l3],'Predicted Matches - Recovered Params',...
            'Predicted Matches - Simulated Params','Simulated Matches');
        lgd.Location = 'northwest';
        xlim([0 1]);
        ylim([0 0.5]);
    end
    % Edit titles if plots were only made for best and worst observers
    if ~p.Results.makeAllObserverPlots
        subplot('Position', subplotPosVectors(2,1).v);
        theTitle = sprintf('L and M Cone Spectral Sensitivities, Observer %g (Best)',bestObs);
        title(theTitle);
        
        subplot('Position', subplotPosVectors(2,2).v);
        theTitle = sprintf('Generalized Pitt Diagram, Observer %g (Best)',bestObs);
        title(theTitle);
        
        subplot('Position', subplotPosVectors(3,1).v);
        theTitle = sprintf('L and M Cone Spectral Sensitivities, Observer %g (Worst)',worstObs);
        title(theTitle);
        
        subplot('Position', subplotPosVectors(3,2).v);
        theTitle = sprintf('Generalized Pitt Diagram, Observer %g (Worst)',worstObs);
        title(theTitle);
    end
    % Save plot
    if isempty(p.Results.plotTitle)
        sgtitle([subjID ' Error'],'Interpreter', 'none');
    else
        sgtitle([p.Results.plotTitle ': Error'],'Interpreter', 'none');
    end
    NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_errPlots.pdf']),...
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
        'topMargin',      0.1);
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
        theTitle = sprintf('%s Recovered vs Simulated',cell2mat(coneParamNames(ii)));
        title(theTitle);
        xlabel('Simulated Parameters');
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
    end
    % Save figure
    if isempty(p.Results.plotTitle)
        sgtitle([subjID ' Parameters'],'Interpreter', 'none');
    else
        sgtitle([p.Results.plotTitle ': Parameters'],'Interpreter', 'none');
    end
    NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_paramPlots.pdf']),...
        paramsFig,300);
end
end