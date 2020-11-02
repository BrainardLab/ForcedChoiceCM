function [questData,psiParamsQuest,psiParamsFit] = ...
    qpRayleighSim(subjID,nObservers,nTrials,baseConeParams,coneParamsToVary,....
    noiseScaleFactor,p1Wl,p2Wl,testWls,varargin)
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
%                         [0.8078 4.1146 1.2592 0.02].
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'lambdaRef'         -Numerical scale factor between 0 and 1 for the
%                         value of lambda used for the reference light,
%                         which is used as a baseline for computing
%                         opponent contrasts. Default is 0.8.
%    'sampledObservers'  -nObservers x 8 array of previously-sampled
%                         observer parameters, useful for when observers
%                         are used across multiple experimental conditions.
%                         Default is [].
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'plotAll'           -Logical. If true, make plots summarizing QUEST+
%                         results for all subjects. Default is false.
%    'plotLast'           -Logical. If true, make plots summarizing QUEST+
%                         results for the last subject. Default is true.

% History
%    10/13/20   dce   -Wrote it
%    10/21/20   dce   -Changed to use a single QUEST object
%    10/23/20   dce   -Added likelihood function and plotting
%    10/26/20   dhb   -Commenting and a little cosmetic editing.
%    10/27/20   dce   -Allowed multiple observers, changed parameters to
%                      function inputs
%    10/28/20   dce   -Modified plotting

close all;

% Parse optional input
p = inputParser;
p.addParameter('plotAll',false,@(x)(islogical(x)));
p.addParameter('plotLast',true,@(x)(islogical(x)));
p.addParameter('precomputeQuest',true,@(x)(islogical(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[0.8078 4.1146 1.2592 0.0200],@(x)(isvector(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('lambdaRef',0.8,@(x)(isnumeric(x)));
p.addParameter('sampledObservers',[],@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Define output directory 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'Quest',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
fName = [subjID '.mat'];
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

% Set up the standard observer, and extract its base cone parameters and
% opponent parameters.
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
psiParamsDomainList = {linspace(-2,2,10),linspace(-2,2,10),...
    linspace(-2,2,10),linspace(-2,2,10),linspace(-2,2,10),...
    linspace(-2,2,9),linspace(-2,2,9),linspace(-2,2,9)};
for i = 1:length(indDiffSds)
    psiParamsDomainList{i} = psiParamsDomainList{i}*coneParamsToVary(i)*indDiffSds(i);
    psiParamsDomainList{i} = unique(psiParamsDomainList{i});
end

% Stimulus parameter list.
% The first parameter is lambda ( = proportion red in the primary mixture),
% the second is test intensity, and the third is the test wavelength.
stimParamsDomainList = {0:0.05:1,0:0.05:1,testWls};

% Psychometric function for Rayleigh matching
PFSim = @(stimParams,coneParams)qpPFRMFull(stimParams,coneParams,...
    p.Results.opponentParams,p.Results.opponentParams(4)*noiseScaleFactor,...
    p.Results.S,p1Spd,p2Spd,testSpds,testWls,p.Results.lambdaRef);

% Set up a Quest object, with an option to use precomputed data if
% available
if (~p.Results.precomputeQuest)
    startTime = tic;
    fprintf('Initializing quest structure ...\n');
    questDataRaw = qpInitialize('nOutcomes', 4, ...
        'qpPF',PFSim,'qpOutcomeF',[],...
        'stimParamsDomainList',stimParamsDomainList, ...
        'psiParamsDomainList',psiParamsDomainList, ...
        'verbose', true);
    elapsedTime = toc(startTime);
    stimParamsDomainListCheck = stimParamsDomainList;
    psiParamsDomainListCheck = psiParamsDomainList;
    fprintf('Done initializing in %0.1f seconds\n',elapsedTime);
    save('questDataRaw','questDataRaw','stimParamsDomainListCheck','psiParamsDomainListCheck','-v7.3');
else
    fprintf('Loading quest structure ...\n');
    load questDataRaw questDataRaw stimParamsDomainListCheck psiParamsDomainListCheck;
    if (length(stimParamsDomainList) ~= length(stimParamsDomainListCheck))
        error('Change in stim parameters since cache of questDataRaw');
    end
    for ii = 1:length(stimParamsDomainList)
        if (any(stimParamsDomainList{ii} ~= stimParamsDomainListCheck{ii}))
            error('Change in stim parameters since cache of questDataRaw');
        end
    end
    if (length(psiParamsDomainList) ~= length(psiParamsDomainListCheck))
        error('Change in psi parameters since cache of questDataRaw');
    end
    for ii = 1:length(psiParamsDomainList)
        if (any(psiParamsDomainList{ii} ~= psiParamsDomainListCheck{ii}))
            error('Change in psi parameters since cache of questDataRaw');
        end
    end
    fprintf('Done loading\n');
end

% Initialize data arrays to store results for each observer
questData = cell(1,nObservers);                            % QUEST objects
psiParamsQuest = zeros(nObservers,length(baseConeParams)); % Recovered params (max posterior)
psiParamsFit = zeros(nObservers,length(baseConeParams));   % Recovered params (maximum likelihood)
trialPrintout = 20;  % Print message after every 20 trials

% Nominal matches for each observer and test wavelength, in the form [lambda testIntensity]
nominalMatch = zeros(nObservers,length(testWls),2); 

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
            fprintf('\t Observer %d of %d, Trial %d of %d\n',...
                ii,nObservers,tt,nTrials);
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
        psiParamsQuest(ii,:),questData{ii}.nOutcomes,'lowerBounds',...
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
    
    % Make plots, if desired
    if p.Results.plotAll || (p.Results.plotLast && ii == nObservers)
        % Plot of QUEST runs
        figure(); clf;
        stimCounts = qpCounts(qpData(questData{ii}.trialData),questData{ii}.nOutcomes);
        lambdaVals = [];
        tIVals = [];
        tWlVals = [];
        markerColors = [];
        markerSizes = [];
        for cc = 1:length(stimCounts)
            stim = stimCounts(cc).stim;
            for jj = 1:questData{ii}.nOutcomes
                outcomeCount = stimCounts(cc).outcomeCounts(jj);
                if outcomeCount == 0
                    continue;
                end
                lambdaVals = [lambdaVals, stim(1)];
                tIVals = [tIVals, stim(2)];
                tWlVals = [tWlVals, stim(3)];
                markerColors = [markerColors, jj];
                markerSizes = [markerSizes,1000*outcomeCount/max(nTrials)];
            end
        end
        scatter3(lambdaVals,tIVals,tWlVals,markerSizes,markerColors,'o','LineWidth',2,...
            'MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4);
        hold on;
        scatter3(nominalMatch(ii,:,1),nominalMatch(ii,:,2),testWls,'k*');
        hold off;
        xlim([0 1]);
        ylim([0 1]);
        zlim([min(testWls)-10,max(testWls)+10]);
        xlabel('Lambda (Proportion Red)');
        ylabel('Test Intensity');
        zlabel('Peak Test Wavelength (nm)');
        title('QUEST+ Trial Placement');
        legend('QUEST+ Trials', 'Nominal Matches'); 
        
        % Plot of cone parameters
        paramsFig = figure(); clf; 
        set(paramsFig,'Color',[1 1 1],'Position',[10 10 1700 800]);
        hold on;
        nCols2 = 4;
        nRows2 = ceil(length(psiParamsFit(ii,:))/nCols2);
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
        for kk = 1:length(psiParamsFit(ii,:))
            % Make a subplot in the correct position
            row = ceil(kk/nCols2);
            col = mod(kk,nCols2);
            if col == 0
                col = nCols2;
            end
            subplot('Position', subplotPosVectors(row,col).v);
            hold on;
        
            % Define axis limits
            if (kk==6) || (kk==7) || (kk==8)  % Lambda max shifts, in nm
                limits = [-5 5];
            else                             % Density shifts, in percent
                limits = [-40 40];
            end
            xlim(limits);
            ylim(limits);
            axis('square');
        
            % Plot data
            lambdaVals = simConeParams(kk);   % Predicted parameters
            tIVals = psiParamsFit(ii,kk); % Recovered params
            plot(lambdaVals,tIVals,'b* ','MarkerSize',7,'LineWidth',1);
            refline(1,0);
        
            % Titles and labels
            theTitle = sprintf('%s Recovered vs Simulated',cell2mat(coneParamNames(kk)));
            title(theTitle);
            xlabel('Simulated Parameters');
            ylabel('Recovered Parameters');
        end
        
        % Plots of cone fundamentals
        observerRecovered = genRayleighObserver('coneVec',psiParamsFit(ii,:),...
            'opponentParams',p.Results.opponentParams,'S',p.Results.S);
        figure(); clf; 
        hold on;
        h1 = plot(wls,observer.T_cones(1:2,:),'r-','LineWidth',2.5);
        h2 = plot(wls,observerRecovered.T_cones(1:2,:),'b-','LineWidth',1.5);
        legend([h1(1) h2(1)],'Simulated Cones', 'Recovered Cones');
        title('L and M Cones');
        xlabel('Wavelength (nm)');
        ylabel('Sensitivity');
        
        figure(); clf;
        hold on;
        subplot(1,2,1); hold on
        plot(wls,observer.T_cones(1,:)-observerRecovered.T_cones(1,:),'r-',...
            'LineWidth',2.5);
        plot(wls,observer.T_cones(1,:)-observerStandard.T_cones(1,:),'k-',...
            'LineWidth',2.5);
        refline(0,0);
        title('L Cone Sensitivity Differences');
        xlabel('Wavelength (nm)');
        ylabel('Sensitivity Difference');
        ylim([-0.06 0.06]);
        legend({'Recovered', 'Standard'},'Location','SouthEast');
        
        subplot(1,2,2); hold on
        plot(wls,observer.T_cones(2,:)-observerRecovered.T_cones(2,:),'r-',...
            'LineWidth',2.5);
        plot(wls,observer.T_cones(2,:)-observerStandard.T_cones(2,:),'k-',...
            'LineWidth',2.5);
        refline(0,0);
        title('M Cone Sensitivity Differences');
        xlabel('Wavelength (nm)');
        ylabel('Sensitivity Difference');
        ylim([-0.06 0.06]);
        legend({'Recovered', 'Standard'},'Location','SouthEast');
    end
    
end
save(outputFile,'subjID','nObservers','nTrials','baseConeParams','coneParamsToVary',....
    'noiseScaleFactor','p1Wl','p2Wl','testWls','p','questData','psiParamsQuest',...
    'psiParamsFit','nominalMatch');
fprintf('Finished observer %g of %g\n', ii, nObservers);
end