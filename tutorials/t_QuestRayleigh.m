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

% Basic apparatus parameters
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
simulatedPsiParamsStruct.coneParams.indDiffParams.lambdaMaxShift = [0 0 0]';
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
stimParamsDomainList = {570:30:660, 540:20:560, 660:20:680, -0.09:0.03:0.09, -0.09:0.03:0.09, 0, -0.09:0.03:0.09, -0.09:0.03:0.09, 0};
psiParamsDomainList = {0, 0, -20:5:20, -20:5:20, 0, -4:1:4, -4:1:4, 0, 0.02};
%psiParamsDomainList = {0, 0, 0, 0, 0, -4:1:4,-4:1:4, 0, 0.02};
psiParamsLabels = {'Lens Density', 'Macular Pigment Density', 'L density', 'M density', 'S density', 'L lambda max', 'M lambda max', 'S lambda max', 'Noise'};

% Marginalize?  Costs considerable time as domain grows.  Empty for not.
marginalizeVector = [];

USE_PRECOMPUTE = false;
if (~USE_PRECOMPUTE)
    fprintf('Initializing quest structure ...\n');
    
    startTime = tic;
    questDataRaw = qpInitialize(...
        'nOutcomes', 2, ...
        'qpPF',qpPFFun, ...
        'stimParamsDomainList',stimParamsDomainList, ...
        'psiParamsDomainList',psiParamsDomainList, ...
        'filterStimParamsDomainFun',@(stimParamsVec) qpFCCMStimDomainCheck(stimParamsVec,stimVecType,stimParamsStruct), ...
        'marginalize', marginalizeVector, ...
        'verbose', true ...
        );
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

%% qpRun estimating the parameters, over and over
rng('default'); rng(3008,'twister');
nParamSets = 20;
nRunsPerParamSet = 2;
nTrials = 256;

[domainVlb,domainVub] = qpGetBoundsFromDomainList(psiParamsDomainList);
for ss = 1:nParamSets
    % Set up parameters for this run by random draw.  But keep noise
    % fixed, as we can probably establish that separately.
    simulatedPsiParamsVecCell{ss} = qpDrawFromDomainList(psiParamsDomainList);
    simulatedPsiParamsVecCell{ss}(9) = simulatedPsiParamsVec(9);
    simulatedPsiParamsVec = simulatedPsiParamsVecCell{ss};
    simulatedPsiParamsStruct = ObserverVecToParams(psiVecType,simulatedPsiParamsVec,simulatedPsiParamsStruct);
    
    % Standard QUEST+ simulated observer.  Ths function definition needs to come after we set
    % stimulus parameter vec and before we call it.
    simulatedObserverFun = @(stimParamsVec) qpSimulatedObserver(stimParamsVec,qpPFFun,simulatedPsiParamsVec);
    
    % Do the runs
    for rr = 1:nRunsPerParamSet
        % Simulate run
        fprintf('*** Simluated run %d of %d for parameters set %d of %d:\n',rr,nRunsPerParamSet,ss,nParamSets);
        startTime = tic;
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
        elapsedTime = toc(startTime);
        fprintf('Done with trial simulation, %0.3f calculation time per trial\n',(elapsedTime)/nTrials);
        
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

%% Save simulation
save('simulatedData','-v7.3');

%% Marginalize posterior to improve estimates of individual parameters.
whichParamsToMarginalize = [1, 2, 3, 4, 5, 8, 9];
% figure;
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        % Marginalize
        psiParamsDomain = questData{ss,rr}.psiParamsDomain;
        posterior = questData{ss,rr}.posterior;
        [marginalPosterior{ss,rr},marginalPsiParamsDomain, marginalPsiParamsLabels] = qpMarginalizePosterior(posterior,psiParamsDomain,whichParamsToMarginalize,psiParamsLabels);
        
        % Find max of marginalized posterior
        marginalPsiParamsIndex = qpListMaxArg(marginalPosterior{ss,rr});
        marginalPsiParamsQuest{ss,rr} = marginalPsiParamsDomain(marginalPsiParamsIndex,:);
        
        % % Little plot
        % clf; hold on;
        % plot3(marginalPsiParamsDomain(:,1),marginalPsiParamsDomain(:,2),marginalPosterior{ss,rr},'ro','MarkerFaceColor','r','MarkerSize',6);
        % plot(simulatedPsiParamsVecCell{ss}(6),simulatedPsiParamsVecCell{ss}(7),'ko','MarkerFaceColor','k','MarkerSize',8);
        % plot(marginalPsiParamsDomain(marginalPsiParamsIndex,1),marginalPsiParamsDomain(marginalPsiParamsIndex,2),'go','MarkerFaceColor','g','MarkerSize',6);
        % xlabel('L lambda max'); ylabel('M lambda max');
        % view([37 85]);
        % drawnow;
        % pause;
    end
end

%% Compute simulated and obtained fundamentals
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        % Get the simuated psi parameters and fundamentals
        simulatedParamsStruct = ObserverVecToParams('basic',simulatedPsiParamsVecCell{ss},psiParamsStruct);
        T_simulated{ss,rr} = ComputeObserverFundamentals(simulatedParamsStruct.coneParams,S);
        
        % Get the fit psi parameters and fundamentals
        fitParamsStruct = ObserverVecToParams('basic',psiParamsFit{ss,rr},psiParamsStruct);
        T_fit{ss,rr} = ComputeObserverFundamentals(fitParamsStruct.coneParams,S);
        
        Lsse(ss,rr) = sqrt(sum( (T_simulated{ss,rr}(1,:)-T_fit{ss,rr}(1,:)).^2 ));
        Msse(ss,rr) = sqrt(sum( (T_simulated{ss,rr}(2,:)-T_fit{ss,rr}(2,:)).^2 ));
        Ssse(ss,rr) = sqrt(sum( (T_simulated{ss,rr}(3,:)-T_fit{ss,rr}(3,:)).^2 ));
    end
end

%% Make plots of simulated versus estimated parameters
theParamIndex = 6;
figure; clf; hold on
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        %plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsQuest{ss,rr}(theParamIndex),'bo','MarkerFaceColor','b','MarkerSize',8);
        %plot(simulatedPsiParamsVecCell{ss}(theParamIndex),marginalPsiParamsQuest{ss,rr}(1),'go','MarkerFaceColor','g','MarkerSize',4);
        %plot(simulatedPsiParamsVecCell{ss}(theParamIndex),marginalPsiParamsQuest1{ss,rr},'kx','MarkerFaceColor','k','MarkerSize',8);
    end
end
xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
axis('square');
xlabel('Simulated Lambda Max Shift');
ylabel('Estimated Max Shift');
title('L Cones');

theParamIndex = 3;
figure; clf; hold on
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsQuest{ss,rr}(theParamIndex),'bo','MarkerFaceColor','b','MarkerSize',8);
    end
end
xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
axis('square');
xlabel('Simulated Density Shift');
ylabel('Estimated Density Shift');
title('L Cones');

theParamIndex1 = 6; theParamIndex2 = 3;
minMarkerSize = 1; maxMarkerSize = 16;
figure; clf; hold on
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        markerSize = minMarkerSize  + (maxMarkerSize - minMarkerSize)*Lsse(ss,rr)/(max([Lsse(:) ; Msse(:)])-min([Lsse(:) ; Msse(:)]));
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex1)-psiParamsFit{ss,rr}(theParamIndex1), ...
             simulatedPsiParamsVecCell{ss}(theParamIndex2)-psiParamsFit{ss,rr}(theParamIndex2),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
        
    end
end
xlim([-4 4]);
ylim([-40 40]);
xlabel('Lambda-max delta error');
ylabel('Density delta error');
title(sprintf('L cones, min sse: %0.4f; max sse %0.4f',min([Lsse(:) ; Msse(:)]),max([Lsse(:) ; Msse(:)])));

theParamIndex = 7;
figure; clf; hold on
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        %plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsQuest{ss,rr}(theParamIndex),'bo','MarkerFaceColor','b','MarkerSize',8);
        %plot(simulatedPsiParamsVecCell{ss}(theParamIndex),marginalPsiParamsQuest{ss,rr}(2),'go','MarkerFaceColor','g','MarkerSize',4);
    end
end
xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
axis('square');
xlabel('Simulated Lambda Max Shift');
ylabel('Estimated Max Shift');
title('M Cones');

theParamIndex = 4;
figure; clf; hold on
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsFit{ss,rr}(theParamIndex),'ro','MarkerFaceColor','r','MarkerSize',8);
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex),psiParamsQuest{ss,rr}(theParamIndex),'bo','MarkerFaceColor','b','MarkerSize',8);
    end
end
xlim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
ylim([domainVlb(theParamIndex) domainVub(theParamIndex)]);
plot([domainVlb(theParamIndex) domainVub(theParamIndex)],[domainVlb(theParamIndex) domainVub(theParamIndex)],'k','LineWidth',1);
axis('square');
xlabel('Simulated Density Shift');
ylabel('Estimated Density Shift');
title('M Cones');

theParamIndex1 = 7; theParamIndex2 = 4;
figure; clf; hold on
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet
        markerSize = minMarkerSize  + (maxMarkerSize - minMarkerSize)*Msse(ss,rr)/(max([Lsse(:) ; Msse(:)])-min([Lsse(:) ; Msse(:)]));
        plot(simulatedPsiParamsVecCell{ss}(theParamIndex1)-psiParamsFit{ss,rr}(theParamIndex1), ...
             simulatedPsiParamsVecCell{ss}(theParamIndex2)-psiParamsFit{ss,rr}(theParamIndex2),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
        
    end
end
xlim([-4 4]);
ylim([-40 40]);
xlabel('Lambda-max delta error');
ylabel('Density delta error');
title(sprintf('M cones, min sse: %0.4f; max sse %0.4f',min([Lsse(:) ; Msse(:)]),max([Lsse(:) ; Msse(:)])));

%% Evaluate some specific runs
fundamentalsFig = figure; clf;
set(gcf,'Position',[50 420 1800 900]);
for ss = 1:nParamSets
    for rr = 1:nRunsPerParamSet  
        figure(fundamentalsFig); clf;
        subplot(1,2,1); hold on
        plot(SToWls(S),TRef(1,:),'k:','LineWidth',2);
        plot(SToWls(S),T_simulated{ss,rr}(1,:),'r','LineWidth',3);
        plot(SToWls(S),T_fit{ss,rr}(1,:),'b','LineWidth',2);
        xlabel('Wavelength (nm)');
        ylabel('Fundamental');
        title({sprintf('L cone sim/fit params, l-max: %0.1f/%0.1f; density %0.1f/%0.1f, sse %0.4f', ...
            simulatedPsiParamsVecCell{ss}(6),psiParamsFit{ss,rr}(6),simulatedPsiParamsVecCell{ss}(3),psiParamsFit{ss,rr}(3),Lsse(ss,rr)) ; ' '});
        legend({'Reference', 'Simulated','Fit'});
        
        subplot(1,2,2); hold on
        plot(SToWls(S),TRef(2,:),'k:','LineWidth',2);
        plot(SToWls(S),T_simulated{ss,rr}(2,:),'r','LineWidth',3);
        plot(SToWls(S),T_fit{ss,rr}(2,:),'b','LineWidth',2);
        xlabel('Wavelength (nm)');
        ylabel('Fundamental');
        title({sprintf('M cone sim/fit params, l-max: %0.1f/%0.1f; density %0.1f/%0.1f, sse %0.4f', ...
            simulatedPsiParamsVecCell{ss}(7),psiParamsFit{ss,rr}(7),simulatedPsiParamsVecCell{ss}(4),psiParamsFit{ss,rr}(4),Msse(ss,rr)) ; ' '});
        legend({'Reference', 'Simulated','Fit'});
        
        % subplot(1,3,3); hold on
        % plot(SToWls(S),TRef(3,:),'k:','LineWidth',2);
        % plot(SToWls(S),T_simulated{ss,rr}(3,:),'r','LineWidth',3);
        % plot(SToWls(S),T_fit{ss,rr}(3,:),'b','LineWidth',2);
        % xlabel('Wavelength (nm)');
        % ylabel('Fundamental');
        % title('S cone');
        % legend({'Reference', 'Simulated','Fit'});
        
        drawnow;
        pause;
        commandwindow;
    end
end

%%
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
%
% questParamsStruct = ObserverVecToParams('basic',psiParamsQuest,psiParamsStruct);
%
%
% T_quest = ComputeObserverFundamentals(questParamsStruct.coneParams,S);
%
%


