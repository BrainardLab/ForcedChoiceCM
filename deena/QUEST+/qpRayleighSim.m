function [questData,psiParamsQuest,psiParamsFit] = ...
    qpRayleighSim(subjID,sessionNum,nObservers,nTrials,baseConeParams,....
    coneParamsToVary,noiseScaleFactor,p1Wl,p2Wl,testWls,varargin)
% Runs Quest+ Rayleigh matching experiments for a series of observers.
%
% Syntax:
%   qpRayleighSim(subjID,nObservers,nTrials,baseConeParams,coneParamsToVary,
%   noiseScaleFactor,p1Wl,p2Wl,testWls)
%
% Description:
%    This function uses the QUEST+ routine to conduct a Rayleigh matching
%    experiment and use it to recover observer individual difference
%    parameters.
%
% Inputs:
%    subjID            -Character vector of subject ID
%    sessionNum        -Numeric session ID.
%    nObservers        -Number of simulated observers to test.
%    baseConeParams    -Eight-element numeric vector of individual
%                       difference parameters: lens pigment density,
%                       macular pigment density, L/M/S photopigment
%                       densities, and L/M/S peak spectral sensitivities
%                       (lambda max). The values are used as means for
%                       observer sampling.
%    coneParamsToVary  -Eight-element numeric vector of ones and zeros
%                       indicating which individual difference parameters
%                       should be varied (the noise parameter is excluded).
%                       Parameters set to 1 will be sampled around their
%                       means, while parameters set to 0 will stay at the
%                       values specified in baseParams.
%    noiseScaleFactor  -Number >=0 which determines which scalar multiple
%                       of the opponent noise SD should be used as the
%                       observer noise SD.
%    p1Wl              -Integer value of desired wavelength for the first
%                       primary.
%    p2Wl              -Integer value of desired wavelength for the second
%                       primary.
%    testWls           -Integer or numeric vector of desired wavelengths
%                       for the test light.
%
% Outputs:
%    questData        -1 x nObservers cell array, where each entry contains
%                      the QUEST structure for that observer
%    psiParamsQuest   -nObservers x 8 matrix, where each row contains the
%                      recovered QUEST cone parameters for a different
%                      observer (max posterior estimates).
%    psiParamsFit     -nObservers x 8 matrix, where each row contains the
%                      maximum likelihood estimates of recovered cone
%                      parameters.
%
% Optional key-value pairs:
%    'precomputeQuest'   -Logical. When true, the program uses precomputed
%                         QUEST data if it is available. Default is true.
%    'fieldSize'         -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'opponentParams'    -4-element vector with opponent contrast
%                         parameters. (1) is the luminance weight, (2) is
%                         the RG weight, (3) is the BY weight, and (4) is
%                         the baseline noise standard deviation. Default is
%                         [40.3908  205.7353   62.9590    1.0000].
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'nStimValues'       -Number of possible stimulus values for lambda and
%                         test intensity, spaced evenly between 0 and 1.
%                         Default is 21.
%    'nPsiValues'        -Number of possible levels for the recovered cone
%                         parameters, spaced evenly between -2 and 2
%                         standard deviations. Default is 8.
%    'lambdaRef'         -Numerical scale factor between 0 and 1 for the
%                         value of lambda used for the reference light,
%                         which is used as a baseline for computing
%                         opponent contrasts. Default is 0.8. If set to [],
%                         no reference light is used, and the opponent
%                         contrast of the test is calculated relative to
%                         the primary mixture.
%    'sampledObservers'  -nObservers x 8 array of previously-sampled
%                         observer parameters, useful for when observers
%                         are used across multiple experimental conditions.
%                         Default is [].
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'plotAll'           -Logical. If true, make plots summarizing QUEST+
%                         results for all subjects. Default is false.
%    'plotLast'          -Logical. If true, make plots summarizing QUEST+
%                         results for the last subject. Default is true.
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on
%                         stimulus parameters. Each row represents a given
%                         test wavelength and the limits which are associated
%                         with it. The columns are arranged as follows:
%                         [test wl, min lambda, max lambda, min test
%                         intensity, max test intensity]. Default is [].
% History
%    10/13/20   dce   -Wrote it
%    10/21/20   dce   -Changed to use a single QUEST object
%    10/23/20   dce   -Added likelihood function and plotting
%    10/26/20   dhb   -Commenting and a little cosmetic editing.
%    10/27/20   dce   -Allowed multiple observers, changed parameters to
%                      function inputs
%    10/28/20   dce   -Modified plotting
%    11/4/20    dce   -Added stimulus spacing as a key value pair, edited
%                      plots
%    11/14/20   dce   -Added reference spd calculation
%    11/15/20   dce   -Added option to filter stimuli
%    12/02/20   dce   -Changed how saved QUEST+ data structs are loaded,
%                      added nPsiValues as a key-value pair
%    12/03/20   dce   -Changed how files are saved

% Examples:
%{
    limMatrix = [570.0000   0    0.1264    0.0399    0.0459;...
      590.0000    0.0456    0.4746    0.0462    0.0716;...
      610.0000    0.2617    0.8120    0.0695    0.1325;...
      630.0000    0.6046    0.9604    0.1619    0.2685;...
      650.0000    0.8688    0.9938    0.5109    0.6458];
    qpRayleighSim('testQPRayleigh',1,1,135,zeros(1,8),[0 0 1 1 0 1 1 0],0,670,560,...
        [570 590 610 630 650],'precomputeQuest',false,'nStimValues',101,'stimLimits',limMatrix);
%}
close all;

% Parse optional input
p = inputParser;
p.addParameter('plotAll',false,@(x)(islogical(x)));
p.addParameter('plotLast',true,@(x)(islogical(x)));
p.addParameter('precomputeQuest',true,@(x)(islogical(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('nStimValues',21,@(x)(isnumeric(x)));
p.addParameter('nPsiValues',8,@(x)(isnumeric(x)));
p.addParameter('lambdaRef',0.8,@(x)(isnumeric(x)));
p.addParameter('sampledObservers',[],@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Define output directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'Quest',subjID,num2str(sessionNum));
if exist(outputDir,'dir')
    error('Subject/session directory already exists');
end 
mkdir(outputDir);
fName = [subjID '_' num2str(sessionNum) '_QUEST.mat'];
outputFile = fullfile(outputDir,fName);

% Define spectra
wls = SToWls(p.Results.S);

p1Spd = zeros(length(wls),1);
p1Spd(wls==p1Wl) = 1*p.Results.p1Scale;

p2Spd = zeros(length(wls),1);
p2Spd(wls==p2Wl) = 1*p.Results.p2Scale;

testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1*p.Results.testScale;
end

% Set up the standard observer
observerStandard = genRayleighObserver('age',p.Results.age,...
    'fieldSize',p.Results.fieldSize,'S',p.Results.S,'opponentParams',...
    p.Results.opponentParams,'coneVec',baseConeParams);

% Get cone parameters for the simulated observers. These are either passed
% in as a key-value pair or are sampled using sampleRayleighObservers
if isempty(p.Results.sampledObservers)
    sampledConeParams = sampleRayleighObservers(nObservers,baseConeParams,...
        coneParamsToVary);
else
    sampledConeParams = p.Results.sampledObservers;
    [r,c] = size(sampledConeParams);
    if r~=nObservers || c~=length(baseConeParams)
        error('Provided observer params have incorrect dimensions');
    end
end

% Spacing for possible individual difference values in terms of stddev multiple.
% Note that lambda maxes are slightly more coarsly spaced than optical densities.
indDiffSds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3]; % Parameter standard deviations
psiParamsDomainList = {linspace(-2,2,p.Results.nPsiValues),linspace(-2,2,p.Results.nPsiValues),...
    linspace(-2,2,p.Results.nPsiValues),linspace(-2,2,p.Results.nPsiValues)...
    ,linspace(-2,2,p.Results.nPsiValues),linspace(-2,2,p.Results.nPsiValues),...
    linspace(-2,2,p.Results.nPsiValues),linspace(-2,2,p.Results.nPsiValues)};
for i = 1:length(indDiffSds)
    psiParamsDomainList{i} = psiParamsDomainList{i}*coneParamsToVary(i)*indDiffSds(i);
    psiParamsDomainList{i} = unique(psiParamsDomainList{i});
end

% Stimulus parameter list.
% The first parameter is lambda ( = proportion red in the primary mixture),
% the second is test intensity, and the third is the test wavelength.
stimParamsDomainList = {linspace(0,1,p.Results.nStimValues),...
    linspace(0,1,p.Results.nStimValues),testWls};

% Stimulis parameter filtering function
if isempty(p.Results.stimLimits)
    stimParamFilterFun = [];
else
    stimParamFilterFun = @(stimParams) qpRayleighLimitParams(stimParams,p.Results.stimLimits);
end

% Define a reference spectrum if available
if isempty(p.Results.lambdaRef)
    refSpd = [];
else
    refSpd = p.Results.lambdaRef*p1Spd + (1-p.Results.lambdaRef)*p2Spd;
end

% Psychometric function for Rayleigh matching
PFSim = @(stimParams,coneParams)qpPFRMFull(stimParams,coneParams,...
    p.Results.opponentParams,p.Results.opponentParams(4)*noiseScaleFactor,...
    p.Results.S,p1Spd,p2Spd,testSpds,testWls,'refSpd',refSpd);
nOutcomes = 4;

% Load precomputed QUEST+ data if available, or generate a new file
questDataRaw = [];
qpDataFName = sprintf('questDataRaw_%g_%g.mat',p.Results.nStimValues,...
    p.Results.nPsiValues);
qpDataDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'Quest','qpDataStructs',qpDataFName);
if p.Results.precomputeQuest
    if exist(qpDataDir,'file')
        fprintf('Loading quest structure ...\n');
        theData = load(qpDataDir);
        
        % Input checking
        for ii = 1:length(stimParamsDomainList)
            if (any(stimParamsDomainList{ii} ~= theData.stimParamsDomainListCheck{ii}))
                disp('Error');
                error('Loaded stimulus parameters do not match current settings');
            end
        end
        for ii = 1:length(psiParamsDomainList)
            if (any(psiParamsDomainList{ii} ~= theData.psiParamsDomainListCheck{ii}))
                disp('Error');
                error('Loaded psi parameters do not match current settings');
            end
        end
        if (theData.p1Wl~=p1Wl || theData.p2Wl~=p2Wl ||...
                theData.p.Results.age~=p.Results.age || ...
                theData.p.Results.fieldSize~=p.Results.fieldSize || ...
                ~all(theData.p.Results.opponentParams==p.Results.opponentParams) || ...
                theData.p.Results.p1Scale~=p.Results.p1Scale || ...
                theData.p.Results.p2Scale~=p.Results.p2Scale || ...
                theData.p.Results.testScale~=p.Results.testScale || ...
                theData.p.Results.lambdaRef~=p.Results.lambdaRef || ...
                ~all(theData.p.Results.stimLimits==p.Results.stimLimits))
            
            disp('Error');
            error('Loaded experimental settings do not match current settings');
        end
        questDataRaw = theData.questDataRaw;
    end
end
if (~p.Results.precomputeQuest || isempty(questDataRaw))
    startTime = tic;
    fprintf('Initializing quest structure ...\n');
    questDataRaw = qpInitialize('nOutcomes', nOutcomes, ...
        'qpPF',PFSim,'qpOutcomeF',[],...
        'stimParamsDomainList',stimParamsDomainList, ...
        'psiParamsDomainList',psiParamsDomainList, ...
        'filterStimParamsDomainFun',stimParamFilterFun,'verbose', true);
    elapsedTime = toc(startTime);
    stimParamsDomainListCheck = stimParamsDomainList;
    psiParamsDomainListCheck = psiParamsDomainList;
    fprintf('Done initializing in %0.1f seconds\n',elapsedTime);
    save(qpDataDir,'questDataRaw','stimParamsDomainListCheck',...
        'psiParamsDomainListCheck','p','p1Wl','p2Wl','-v7.3');
end

% Initialize data arrays to store results for each observer
questData = cell(1,nObservers);                            % QUEST objects
psiParamsQuest = zeros(nObservers,length(baseConeParams)); % Recovered params (max posterior)
psiParamsFit = zeros(nObservers,length(baseConeParams));   % Recovered params (maximum likelihood)

% Nominal matches for each observer and test wavelength, in the form [lambda testIntensity]
nominalMatch = zeros(nObservers,length(testWls),2);
trialData = cell(1,nObservers);
trialPrintout = 20;  % Print message after every 20 trials

% Loop through observers
for ii = 1:nObservers
    % Generate simulated observer and associated function
    simConeParams = sampledConeParams(ii,:); % Cone vector for the observer
    % Check that we are not varying a locked parameter
    for kk = 1:length(simConeParams)
        if (simConeParams(kk) ~= 0 && coneParamsToVary(kk) == 0)
            error('Varying a simulated parameter that we are telling Quest+ is locked');
        end
    end
    observer = genRayleighObserver('age',p.Results.age,...
        'fieldSize',p.Results.fieldSize,'S',p.Results.S,'opponentParams',...
        p.Results.opponentParams,'coneVec',simConeParams);
    simObserverFun = @(stimParams) qpSimulatedObserver(stimParams,PFSim,simConeParams);
    
    % Initialize Quest data for this run
    questData{ii} = questDataRaw;
    
    % Run trials
    startTime = tic;
    for tt = 1:nTrials
        % Get stimulus for this trial
        stim = qpQuery(questData{ii});
        
        % Simulate outcome
        outcome = simObserverFun(stim);
        
        % Update quest data structure
        questData{ii} = qpUpdate(questData{ii},stim,outcome);
        
        % Periodic printout
        if (rem(tt,trialPrintout) == 0)
            fprintf('\t Observer %d of %d, Trial %d of %d, time = %d\n',...
                ii,nObservers,tt,nTrials,toc(startTime));
        end
    end
    elapsedTime = toc(startTime);
    fprintf('Done with observer %d simulation, %0.3f calculation time per trial\n',...
        ii,(elapsedTime)/nTrials);
    
    % Max posterior estimate
    psiParamsIndex = qpListMaxArg(questData{ii}.posterior);
    psiParamsQuest(ii,:) = questData{ii}.psiParamsDomain(psiParamsIndex,:);
    fprintf('Simulated cone parameters: %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f \n', ...
        simConeParams(1),simConeParams(2),simConeParams(3),simConeParams(4), ...
        simConeParams(5),simConeParams(6),simConeParams(7),simConeParams(8));
    fprintf('Max posterior QUEST+ cone parameters: %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f \n', ...
        psiParamsQuest(ii,1),psiParamsQuest(ii,2),psiParamsQuest(ii,3),...
        psiParamsQuest(ii,4),psiParamsQuest(ii,5),psiParamsQuest(ii,6),...
        psiParamsQuest(ii,7),psiParamsQuest(ii,8));
    
    % Maximum likelihood fit.
    % Use psiParams from QUEST+ as the starting parameter for the search, and
    % impose as parameter bounds the range provided to QUEST+.
    [domainVlb,domainVub] = qpGetBoundsFromDomainList(psiParamsDomainList);
    psiParamsFit(ii,:) = qpFit(questData{ii}.trialData,questData{ii}.qpPF,...
        psiParamsQuest(ii,:),nOutcomes,'lowerBounds',...
        domainVlb,'upperBounds',domainVub);
    fprintf('Maximum likelihood fit parameters:   %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f\n', ...
        psiParamsFit(ii,1),psiParamsFit(ii,2),psiParamsFit(ii,3),...
        psiParamsFit(ii,4),psiParamsFit(ii,5),psiParamsFit(ii,6),...
        psiParamsFit(ii,7),psiParamsFit(ii,8));
    
    % Compute nominal matches for the observer analytically
    for ww = 1:length(testWls)
        [~,~,nominalMatch(ii,ww,2),nominalMatch(ii,ww,1)] = ...
            computePredictedRayleighMatch(p1Spd,p2Spd,testSpds(:,ww),...
            observer,'addDarkSpd',false);
    end
    trialData{ii} = questData{ii}.trialData;  % Extract trial data to save
    fprintf('Finished observer %g of %g\n', ii, nObservers);
end
save(outputFile,'subjID','nObservers','nTrials','baseConeParams','coneParamsToVary',....
    'noiseScaleFactor','p1Wl','p2Wl','testWls','p','trialData','psiParamsQuest',...
    'psiParamsFit','sampledConeParams','nominalMatch','nOutcomes',...
    'outputDir','sessionNum','-v7.3');

%% Make plots, if desired
if p.Results.plotAll || p.Results.plotLast
    % Plot of cone parameter recovery (all observers). This needs to be
    % edited if we change which parameters we are varying.
    paramsPlot = figure(); clf;
    set(paramsPlot,'Color',[1 1 1],'Position',[10 10 1700 800]);
    hold on;
    nCols2 = 2;
    nRows2 = 2;
    coneParamNamesFull = {'Lens Density','Macular Pigment Density',...
        'L Photopigment Density','M Photopigment Density',...
        'S photopigment density','L Lambda Max','M Lambda Max', 'S Lambda Max'};
    coneParamNames = {'L Photopigment Density','M Photopigment Density',...
        'L Lambda Max','M Lambda Max'};
    paramIndsToPlot = [3 4 6 7];
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', nRows2, ...
        'colsNum', nCols2, ...
        'heightMargin',  0.07, ...
        'widthMargin',    0.07, ...
        'leftMargin',     0.04, ...
        'rightMargin',    0.04, ...
        'bottomMargin',   0.07, ...
        'topMargin',      0.1);

    for kk = 1:nRows2*nCols2
        % Make a subplot in the correct position
        row = ceil(kk/nCols2);
        col = mod(kk,nCols2);
        if col == 0
            col = nCols2;
        end
        subplot('Position', subplotPosVectors(row,col).v);
        hold on;
        
        % Define axis limits
        if (kk==3) || (kk==4)  % Lambda max shifts, in nm
            limits = [-5 5];
        else                   % Density shifts, in percent
            limits = [-40 40];
        end
        xlim(limits);
        ylim(limits);
        axis('square');
        
        % Plot data
        xVals = sampledConeParams(:,paramIndsToPlot(kk));   % Predicted parameters
        yVals = psiParamsFit(:,paramIndsToPlot(kk)); % Recovered params
        plot(xVals,yVals,'b* ','MarkerSize',7);
        refline(1,0);
        
        % Titles and labels
        theTitle = sprintf('%s Recovered vs Simulated',cell2mat(coneParamNames(kk)));
        title(theTitle);
        xlabel('Simulated Parameters');
        ylabel('Recovered Parameters');
    end
    NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_'....
        num2str(sessionNum) '_paramsPlot.pdf']),...
        paramsPlot,300);
  
    % Select observers to plot individual results for
    if p.Results.plotAll
        plotInds = 1:nObservers;
    else
        plotInds = nObservers; % Plot only the last observer
    end
    observerStandard = genRayleighObserver('age',p.Results.age,...
    'fieldSize',p.Results.fieldSize,'S',p.Results.S,'opponentParams',...
    p.Results.opponentParams,'coneVec',baseConeParams);
    for pp = plotInds
        % Plot of QUEST runs for selected observers
        qRunPlot = figure(); clf;
        stimCounts = qpCounts(qpData(trialData{pp}),nOutcomes);
        xVals = zeros(1,length(stimCounts)*nOutcomes);
        yVals = zeros(1,length(stimCounts)*nOutcomes);
        tWlVals = zeros(1,length(stimCounts)*nOutcomes);
        markerColors = zeros(1,length(stimCounts)*nOutcomes);
        markerSizes = zeros(1,length(stimCounts)*nOutcomes);
        arrInd = 0;
        for cc = 1:length(stimCounts)
            stim = stimCounts(cc).stim;
            for jj = 1:nOutcomes
                outcomeCount = stimCounts(cc).outcomeCounts(jj);
                if outcomeCount == 0
                    continue;
                end
                arrInd = arrInd + 1;
                xVals(arrInd) = stim(1);
                yVals(arrInd) = stim(2);
                tWlVals(arrInd) = stim(3);
                markerColors(arrInd) = jj;
                markerSizes(arrInd) = 1000*outcomeCount/max(nTrials);
            end
        end
        % Remove extra zeros from arrays
        xVals = xVals(1:arrInd);
        yVals = yVals(1:arrInd);
        tWlVals = tWlVals(1:arrInd);
        markerColors = markerColors(1:arrInd);
        markerSizes = markerSizes(1:arrInd);
        
        scatter3(xVals,yVals,tWlVals,markerSizes,markerColors,'o','LineWidth',2,...
            'MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4);
        hold on;
        scatter3(nominalMatch(pp,:,1),nominalMatch(pp,:,2),testWls,'k*');
        hold off;
        xlim([0 1]);
        ylim([0 1]);
        zlim([min(testWls)-10,max(testWls)+10]);
        xlabel('Lambda (Proportion Red)');
        ylabel('Test Intensity');
        zlabel('Peak Test Wavelength (nm)');
        title('QUEST+ Trial Placement');
        legend('QUEST+ Trials', 'Nominal Matches');
        NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_' num2str(pp) '_qPlot.pdf']),...
            qRunPlot,300);
        
        % Plot of cone fundamentals
        observerSim = genRayleighObserver('age',p.Results.age,...
            'fieldSize',p.Results.fieldSize,'S',p.Results.S,'opponentParams',...
            p.Results.opponentParams,'coneVec',sampledConeParams(pp,:));
        observerRecovered = genRayleighObserver('coneVec',psiParamsFit(pp,:),...
            'opponentParams',p.Results.opponentParams,'S',p.Results.S,...
            'age',p.Results.age,'fieldSize',p.Results.fieldSize);
        conePlot = figure(); clf;
        wls = SToWls(p.Results.S);
        hold on;
        subplot(3,1,1); hold on
        h1 = plot(wls,observerSim.T_cones(1:2,:),'r-','LineWidth',2.5);
        h2 = plot(wls,observerRecovered.T_cones(1:2,:),'b-','LineWidth',1.5);
        legend([h1(1) h2(1)],'Simulated Cones', 'Recovered Cones');
        title('L and M Cones');
        xlabel('Wavelength (nm)');
        ylabel('Sensitivity');
        
        subplot(3,1,2); hold on
        plot(wls,observerSim.T_cones(1,:)-observerRecovered.T_cones(1,:),'r-',...
            'LineWidth',2.5);
        plot(wls,observerSim.T_cones(1,:)-observerStandard.T_cones(1,:),'k-',...
            'LineWidth',2.5);
        refline(0,0);
        title('L Cone Sensitivity Differences');
        xlabel('Wavelength (nm)');
        ylabel('Sensitivity Difference');
        ylim([-0.1 0.1]);
        legend({'Recovered', 'Standard'},'Location','SouthEast');
        
        subplot(3,1,3); hold on
        plot(wls,observerSim.T_cones(2,:)-observerRecovered.T_cones(2,:),'r-',...
            'LineWidth',2.5);
        plot(wls,observerSim.T_cones(2,:)-observerStandard.T_cones(2,:),'k-',...
            'LineWidth',2.5);
        refline(0,0);
        title('M Cone Sensitivity Differences');
        xlabel('Wavelength (nm)');
        ylabel('Sensitivity Difference');
        ylim([-0.1 0.1]);
        legend({'Recovered', 'Standard'},'Location','SouthEast');
        NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_' num2str(pp) '_conePlot.pdf']),...
            conePlot,300);
    end
end
end