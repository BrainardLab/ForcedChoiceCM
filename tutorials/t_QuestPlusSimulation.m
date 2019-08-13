function t_QuestPlusSimulation
% Demonstrate/test QUEST+ at work on parametric forced choice color matching
%
% Description:
%    This script shows QUEST+ employed to estimate observer parameters for
%    cone fundamentals, based on forced choice judgments.
%
% See also:
%

% History:
%   08/12/19  dhb  Created.

%% Close out stray figures
close all;

%% Set up parameters
%
% Wavelength sampling
S = [400 1 301];

% Apparatus parameters
stimParamsStruct.matchApparatusParams = DefaultMatchApparatusParams('monochromatic',S);
stimParamsStruct.testParams = DefaultTestParams('monochromatic',S);

% Observer parameters
psiParamsStruct.coneParams = DefaultConeParams('cie_asano');
psiParamsStruct.colorDiffParams = DefaultColorDiffParams('opponentContrast');
psiParamsStructRef = psiParamsStruct;

%% Set an adaptation spectrum
adaptationSpd = stimParamsStruct.matchApparatusParams.primaryBasis*[0.25 0.25 0.25]';

%% Define psychometric function in terms of lookup table
stimVectorType = 'basic';
psiVectorType = 'basic';

%% Set up simulated observer parameters
simulatedPsiParamsStruct = psiParamsStruct;
simulatedPsiParamsStruct.coneParams.indDiffParams.dlens = 0;
simulatedPsiParamsStruct.coneParams.indDiffParams.dmac = 0;
simulatedPsiParamsStruct.coneParams.indDiffParams.dphotopigment = [0 0 0]';
simulatedPsiParamsStruct.coneParams.indDiffParams.lambdaMaxShift = [-3 0 0]';
simulatedPsiParamsStruct.colorDiffParams.noiseSd = 0.02;
simulatedPsiParamsVec = ObserverParamsToVec(psiVectorType,simulatedPsiParamsStruct);

%% Psychometric and simulated observer functions
qpPFFun = @(stimParamsVec,psiParamsVec) qpPFFCCM(stimParamsVec,psiParamsVec,S,stimVectorType,stimParamsStruct,psiVectorType,psiParamsStruct,psiParamsStructRef,adaptationSpd);
simulatedObserverFun = @(stimParamsVec) qpSimulatedObserver(stimParamsVec,qpPFFun,simulatedPsiParamsVec);

% Initialize quest data
fprintf('Initializing quest structure ...');
questData = qpInitialize(...
    'nOutcomes', 2, ...
    'qpPF',qpPFFun, ...
    'stimParamsDomainList',{[440 550 660], -0.04:0.04:0.04, 0, 0, 0, -0.04:0.04:0.04, 0}, ...
    'psiParamsDomainList',{0, 0, 0, 0, 0, -4:4, 0, 0, 0.02} ...
    );
% 'filterPsiParamsDomainFun',@(psiParams) qpQuestPlusColorMaterialCubicModelParamsCheck(psiParams,maxStimValue,maxPosition,minSpacing) ...
fprintf(' done\n');      

%% qpRun estimating the parameters
fprintf('*** Simluated run, estimate parametric cone fundamentals:\n');
rng('default'); rng(3008,'twister');
nTrials = 256;
for tt = 1:nTrials
    % Get stimulus for this trial
    stim = qpQuery(questData);
    
    % Simulate outcome
    outcome = simulatedObserverFun(stim);
    
    % Update quest data structure
    questData = qpUpdate(questData,stim,outcome); 
    
    if (rem(tt,10) == 0)
        fprintf('\tTrial %d of %d\n',tt,nTrials);
    end
end
fprintf('Done with trial simulation\n');

%% Process simulated data
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Simulated parameters:            %0.1f, %0.1f, %0.1f, %0.1f. %0.1f, %0.1f, %0.1f, %0.1f, %0.3f\n', ...
    simulatedPsiParamsVec(1),simulatedPsiParamsVec(2),simulatedPsiParamsVec(3),simulatedPsiParamsVec(4), ...
    simulatedPsiParamsVec(5),simulatedPsiParamsVec(6),simulatedPsiParamsVec(7),simulatedPsiParamsVec(8),simulatedPsiParamsVec(9));
fprintf('Max posterior QUEST+ parameters: %0.1f, %0.1f, %0.1f, %0.1f. %0.1f, %0.1f, %0.1f, %0.1f, %0.3f\n', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4), ...
    psiParamsQuest(5),psiParamsQuest(6),psiParamsQuest(7),psiParamsQuest(8),psiParamsQuest(9));
% psiParamsCheck = [-360000 -500000 12000 0];
% assert(all(psiParamsCheck == round(10000*psiParamsQuest)),'No longer get same QUEST+ estimate for this case');

%% Maximum likelihood fit.
% Use psiParams from QUEST+ as the starting parameter for the search, and
% impose as parameter bounds the range provided to QUEST+.
psiParamsFit = qpFit(questData.trialData,questData.qpPF,psiParamsQuest,questData.nOutcomes,...
     'lowerBounds', [0, 0, 0, 0, 0, -4, 0, 0, 0.02],'upperBounds',[0, 0, 0, 0, 0, 4, 0, 0, 0.02]);
fprintf('Maximum likelihood fit parameters:  %0.1f, %0.1f, %0.1f, %0.1f. %0.1f, %0.1f, %0.1f, %0.1f, %0.3f\n', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4), ...
    psiParamsFit(5),psiParamsFit(6),psiParamsFit(7),psiParamsFit(8),psiParamsFit(9));
% psiParamsCheck = [-359123 -485262 11325 0];
% assert(all(psiParamsCheck == round(10000*psiParamsFit)),'No longer get same ML estimate for this case');
 
% Plot trial locations together with maximum likelihood fit.
%
% Point transparancy visualizes number of trials (more opaque -> more
% trials), while point color visualizes percent correct (more blue -> more
% correct).
% figure; clf; hold on
% stimCounts = qpCounts(qpData(questData.trialData),questData.nOutcomes);
% stim = zeros(length(stimCounts),questData.nStimParams);
% for cc = 1:length(stimCounts)
%     stim(cc,:) = stimCounts(cc).stim;
%     nTrials(cc) = sum(stimCounts(cc).outcomeCounts);
%     pCorrect(cc) = stimCounts(cc).outcomeCounts(2)/nTrials(cc);
% end
% for cc = 1:length(stimCounts)
%     h = scatter(stim(cc,1),stim(cc,3),100,'o','MarkerEdgeColor',[1-pCorrect(cc) 0 pCorrect(cc)],'MarkerFaceColor',[1-pCorrect(cc) 0 pCorrect(cc)],...
%         'MarkerFaceAlpha',nTrials(cc)/max(nTrials),'MarkerEdgeAlpha',nTrials(cc)/max(nTrials));
% end
% plotSfs = (0:2:40)';
% [~,plotFitThresholds] = qpPFSTCSF(...
%     [plotSfs zeros(size(plotSfs)) zeros(size(plotSfs))], ...
%     psiParamsFit);
% plot(plotSfs,plotFitThresholds,'-','Color',[1 0.2 0.0],'LineWidth',3);
% xlabel('Spatial Frequency (c/deg)');
% ylabel('Contrast (dB)');
% xlim([0 40]); ylim([-50 0]);
% title({'Estimate spatial CSF', ''});
% drawnow;

