function t_QuestRayleigh
% Demonstrate/test QUEST+ at work on parametric forced choice Rayleigh matching
%
% Description:
%    This script shows QUEST+ employed to estimate observer parameters for
%    cone fundamentals, based on forced choice judgments.
%
% See also:
%

% History:
%   08/24/19  dhb  Created.

%% Close out stray figures
close all;

%% Set up parameters
%
% Wavelength sampling
S = [400 1 301];

% Apparatus parameters
stimParamsStruct.matchApparatusParams = DefaultMatchApparatusParams('rayleigh',S);
stimParamsStruct.testParams = DefaultTestParams('rayleigh',S);

% Observer parameters
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
psiParamsStruct.colorDiffParams = DefaultColorDiffParams('opponentContrast');
psiParamsStructRef = psiParamsStruct;
referencePsiParamsVec = ObserverParamsToVec('basic',psiParamsStructRef);

%% Set an adaptation spectrum
adaptationSpd = stimParamsStruct.matchApparatusParams.primaryBasis*[0.25 0.25]';

%% Define psychometric function in terms of lookup table
stimVecType = 'basic';
psiVecType = 'basic';

%% Set up simulated observer parameters
simulatedPsiParamsStruct = psiParamsStruct;
simulatedPsiParamsStruct.coneParams.indDiffParams.dlens = 0;
simulatedPsiParamsStruct.coneParams.indDiffParams.dmac = 0;
simulatedPsiParamsStruct.coneParams.indDiffParams.dphotopigment = [0 0 0]';
simulatedPsiParamsStruct.coneParams.indDiffParams.lambdaMaxShift = [-2.6 2.3 0]';
simulatedPsiParamsStruct.colorDiffParams.noiseSd = 0.02;
simulatedPsiParamsVec = ObserverParamsToVec(psiVecType,simulatedPsiParamsStruct);

%% Psychometric and simulated observer functions
%
% Call into qpPFFCM
TRef = ComputeObserverFundamentals(psiParamsStructRef.coneParams,S);
adaptationLMSRef = TRef*adaptationSpd;
qpPFFun = @(stimParamsVec,psiParamsVec) qpPFFCCM(stimParamsVec,psiParamsVec,S,stimVecType,stimParamsStruct,psiVecType,psiParamsStruct,psiParamsStructRef,adaptationSpd,TRef,adaptationLMSRef);

% If you want to test with a PF that runs fast, uncomment this line.  The
% resulting behavior will be meaningless, however.
% qpPFFun = @(stimParamsVec,psiParamsVec) 0.5;

%% Initialize quest data
%
% From Asano et al. (2016), PlosOne, Table 1.
%   Lens density standard deviation in percent: ~20.
%   Macular density standard deviation in percent: ~40
%   Photopigment density standard deviation in percent: ~15
%   L lambda max standard deviation in nm: ~3
%   M lambda max standard deviation in nm: ~2
%   S lambda max standard deviation in nm: ~1.5
% These are pretty ballpark numbers.

% This can be slow, so best is to precompute and then load on individual
% runs.
stimParamsDomainList = {540:10:670, -0.09:0.03:0.09, -0.09:0.03:0.09, 0, -0.09:0.03:0.09, -0.09:0.03:0.09, 0};
psiParamsDomainList = {0, 0, -20:5:20, -20:5:20, 0, -4:1:4, -4:1:4, 0, 0.02};
psiParamsDomainList = {0, 0, 0, 0, 0, -4:1:4, 0, 0, 0.02};
USE_PRECOMPUTE = false;
if (~USE_PRECOMPUTE)
    fprintf('Initializing quest structure ...\n');
    
    startTime = GetSecs;
    questDataRaw = qpInitialize(...
        'nOutcomes', 2, ...
        'qpPF',qpPFFun, ...
        'stimParamsDomainList',stimParamsDomainList, ...
        'psiParamsDomainList',psiParamsDomainList, ...
        'filterStimParamsDomainFun',@(stimParamsVec) qpFCCMStimDomainCheck(stimParamsVec,stimVecType,stimParamsStruct), ...
        'verbose', true ...
        );
    stopTime = GetSecs;
    stimParamsDomainListCheck = stimParamsDomainList;
    psiParamsDomainListCheck = psiParamsDomainList;
    fprintf('Done initializing in %0.1f seconds\n',stopTime-startTime);
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

%% qpRun estimating the parameters, over and over
rng('default'); rng(3008,'twister');
nParamSets = 10;
nRunsPerParamSet = 5;
nTrials = 1024;

[domainVlb,domainVub] = GetBoundsFromDomainList(psiParamsDomainList);
for ss = 1:nParamSets
    % Set up parameters for this run by random draw.  But keep noise
    % fixed, as we can probably establish that separately.
    simulatedPsiParamsVecCell{ss} = DrawFromDomainList(psiParamsDomainList);
    simulatedPsiParamsVecCell{ss}(9) = simulatedPsiParamsVec(9);
    simulatedPsiParamsVec = simulatedPsiParamsVecCell{ss};
    simulatedPsiParamsStruct = ObserverVecToParams(psiVecType,simulatedPsiParamsVec,simulatedPsiParamsStruct);
    
    % Standard QUEST+ simulated observer
    simulatedObserverFun = @(stimParamsVec) qpSimulatedObserver(stimParamsVec,qpPFFun,simulatedPsiParamsVec);
    
    % Do the runs
    for rr = 1:nRunsPerParamSet
        % Simulate run
        fprintf('*** Simluated run %d of %d for parameters set %d of %d:\n',rr,nRunsPerParamSet,ss,nParamSets);
        startTime = GetSecs;
        questData{ss,rr} = questDataRaw;
        for tt = 1:nTrials
            % Get stimulus for this trial
            stim = qpQuery(questData{ss,rr});
            
            % Simulate outcome
            outcome = simulatedObserverFun(stim);
            
            % Update quest data structure
            questData{ss,rr} = qpUpdate(questData{ss,rr},stim,outcome);
            
            if (rem(tt,128) == 0)
                fprintf('\tTrial %d of %d\n',tt,nTrials);
            end
        end
        stopTime = GetSecs;
        fprintf('Done with trial simulation, %0.3f calculation time per trial\n',(stopTime-startTime)/nTrials);
        
        % Process simulated data
        psiParamsIndex = qpListMaxArg(questData{ss,rr}.posterior);
        psiParamsQuest{ss,rr} = questData{ss,rr}.psiParamsDomain(psiParamsIndex,:);
        fprintf('Simulated parameters:                %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n', ...
            simulatedPsiParamsVec(1),simulatedPsiParamsVec(2),simulatedPsiParamsVec(3),simulatedPsiParamsVec(4), ...
            simulatedPsiParamsVec(5),simulatedPsiParamsVec(6),simulatedPsiParamsVec(7),simulatedPsiParamsVec(8),simulatedPsiParamsVec(9));
        fprintf('Max posterior QUEST+ parameters:     %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n', ...
            psiParamsQuest{ss,rr}(1),psiParamsQuest{ss,rr}(2),psiParamsQuest{ss,rr}(3),psiParamsQuest{ss,rr}(4), ...
            psiParamsQuest{ss,rr}(5),psiParamsQuest{ss,rr}(6),psiParamsQuest{ss,rr}(7),psiParamsQuest{ss,rr}(8),psiParamsQuest{ss,rr}(9));
        
        %% Maximum likelihood fit.
        % Use psiParams from QUEST+ as the starting parameter for the search, and
        % impose as parameter bounds the range provided to QUEST+.
        psiParamsFit{ss,rr} = qpFit(questData{ss,rr}.trialData,questData{ss,rr}.qpPF,psiParamsQuest{ss,rr},questData{ss,rr}.nOutcomes,...
            'lowerBounds', domainVlb,'upperBounds',domainVub);
        fprintf('Maximum likelihood fit parameters:   %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f, %0.3f\n', ...
            psiParamsFit{ss,rr}(1),psiParamsFit{ss,rr}(2),psiParamsFit{ss,rr}(3),psiParamsFit{ss,rr}(4), ...
            psiParamsFit{ss,rr}(5),psiParamsFit{ss,rr}(6),psiParamsFit{ss,rr}(7),psiParamsFit{ss,rr}(8),psiParamsFit{ss,rr}(9));
    end
end

%% Make plots of simulated versus estimated parameters
theParamIndex = 6;
theParamName = 'L Lambda Max Shift';
if (domainVlb(theParamIndex) < domainVub(theParamIndex))
    figure; clf; hold on
    for ss = 1:nParamSets
        for rr = 1:nRunsPerParamSet
            plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        end
    end
    xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
    axis('square');
    xlabel('Simulated');
    ylabel('Estimated');
    title(theParamName);
end

theParamIndex = 7;
theParamName = 'M Lambda Max Shift';
if (domainVlb(theParamIndex) < domainVub(theParamIndex))
    figure; clf; hold on
    for ss = 1:nParamSets
        for rr = 1:nRunsPerParamSet
            plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        end
    end
    xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
    axis('square');
    xlabel('Simulated');
    ylabel('Estimated');
    title(theParamName);
end

theParamIndex = 3;
theParamName = 'L Density';
if (domainVlb(theParamIndex) < domainVub(theParamIndex))
    figure; clf; hold on
    for ss = 1:nParamSets
        for rr = 1:nRunsPerParamSet
            plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        end
    end
    xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
    axis('square');
    xlabel('Simulated');
    ylabel('Estimated');
    title(theParamName);
end

theParamIndex = 4;
theParamName = 'M Density';
if (domainVlb(theParamIndex) < domainVub(theParamIndex))
    figure; clf; hold on
    for ss = 1:nParamSets
        for rr = 1:nRunsPerParamSet
            plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        end
    end
    xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
    plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
    axis('square');
    xlabel('Simulated');
    ylabel('Estimated');
    title(theParamName);
end

%% Evaluate fits a little bit
%
% Get stimulus counts
% stimCounts = qpCounts(qpData(questData.trialData),questData.nOutcomes);
% 
% % Log likelihoods
% logLikelihoodRef = -qpFitError(referencePsiParamsVec,stimCounts,questData.qpPF);
% logLikelihoodSimulated = -qpFitError(simulatedPsiParamsVec,stimCounts,questData.qpPF);
% logLikelihoodQuest = -qpFitError(psiParamsQuest,stimCounts,questData.qpPF);
% logLikelihoodFit = -qpFitError(psiParamsFit,stimCounts,questData.qpPF);
% fprintf('Log likelihoods:\n');
% fprintf('\tReference fundamentals: %0.2f\n',logLikelihoodRef);
% fprintf('\tQuest max posterior:    %0.2f\n',logLikelihoodQuest);
% fprintf('\tSimulated fundamentas:  %0.2f\n',logLikelihoodSimulated);
% fprintf('\tMax likelihood fit:     %0.2f\n',logLikelihoodFit);

% % Compute fundamentals
% simulatedParamsStruct = ObserverVecToParams('basic',simulatedPsiParamsVec,psiParamsStruct);
% questParamsStruct = ObserverVecToParams('basic',psiParamsQuest,psiParamsStruct);
% fitParamsStruct = ObserverVecToParams('basic',psiParamsFit,psiParamsStruct);
% T_simulated = ComputeObserverFundamentals(simulatedParamsStruct.coneParams,S);
% T_quest = ComputeObserverFundamentals(questParamsStruct.coneParams,S);
% T_fit = ComputeObserverFundamentals(fitParamsStruct.coneParams,S);
% 
% fundamentalsFig = figure; clf;
% set(gcf,'Position',[50 420 2400 900]);
% subplot(1,3,1); hold on
% plot(SToWls(S),TRef(1,:),'k:','LineWidth',2);
% plot(SToWls(S),T_simulated(1,:),'r','LineWidth',3);
% %plot(SToWls(S),T_quest(1,:),'b:','LineWidth',2);
% plot(SToWls(S),T_fit(1,:),'b','LineWidth',2);
% xlabel('Wavelength (nm)');
% ylabel('Fundamental');
% title('L cone');
% legend({'Reference', 'Simulated','Fit'});
% 
% subplot(1,3,2); hold on
% plot(SToWls(S),TRef(2,:),'k','LineWidth',3);
% plot(SToWls(S),T_simulated(2,:),'r','LineWidth',3);
% plot(SToWls(S),T_quest(2,:),'b:','LineWidth',2);
% plot(SToWls(S),T_fit(2,:),'b','LineWidth',2);
% xlabel('Wavelength (nm)');
% ylabel('Fundamental');
% title('M cone');
% legend({'Reference', 'Simulated','Fit'});
% 
% subplot(1,3,3); hold on
% plot(SToWls(S),TRef(3,:),'k','LineWidth',3);
% plot(SToWls(S),T_simulated(3,:),'r','LineWidth',3);
% plot(SToWls(S),T_quest(3,:),'b:','LineWidth',2);
% plot(SToWls(S),T_fit(3,:),'b','LineWidth',2);
% xlabel('Wavelength (nm)');
% ylabel('Fundamental');
% title('S cone');
% legend({'Reference', 'Simulated','Fit'});

