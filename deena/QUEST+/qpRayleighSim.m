function qpRayleighSim()
% QUEST+ Rayleigh matching experiment
% Description:
%    This function uses the QUEST+ routine to conduct a Rayleigh matching
%    experiment and use it to recover observer individual difference
%    parameters.

% History
%    10/13/20   dce   -Wrote it
%    10/21/20   dce   -Changed to use a single QUEST object
%    10/23/20   dce   -Added likelihood function and plotting
%    10/26/20   dhb   -Commenting and a little cosmetic editing.

close all; 

% Spectra
% Define monochromatic spectra to use throughout the script
S = [400 1 301];
wls = SToWls(S);

p1Spd = zeros(length(wls),1);
p1Wl = 670;
p1Spd(wls==p1Wl) = 1;

p2Spd = zeros(length(wls),1);
p2Wl = 560;
p2Spd(wls==p2Wl) = 1;

testWls = [570 590 610 630 650]; % Possible wavelengths for the test light
testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1;
end

% Set up the simulated observer.
% coneVec is a vector of eight individual difference parameters:
% lens pigment density, macular pigment density, L/M/S photopigment
% densities, and L/M/S peak spectral sensitivities (lambda max)

% coneVec = sampleRayleighObservers(1,zeros(1,8),[0 0 0 0 0 1 0 0]);
% simConeParams = [0 0 20 10 0 0 0 0];
simConeParams = [0 0 20 10 0 -3 2 0];
observer = genRayleighObserver('coneVec',simConeParams,'S',S);
observerStandard = genRayleighObserver('coneVec',zeros(size(simConeParams)),'S',S);
opponentVec = [observer.colorDiffParams.lumWeight,...
    observer.colorDiffParams.rgWeight,observer.colorDiffParams.byWeight,...
    observer.colorDiffParams.noiseSd];

% Set up psychometric function parameters.  First indicator variable
% determines what parameters we allow to vary.  These are set up in terms
% of a multiple of the population standard deviations.
variableParams = [0 0 1 1 0 1 1 0];               % Which cone params are we varying?
indDiffSds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3]; % Param standard deviations

% Spacing for possible individual difference values in terms of stddev multiple.
% Note that lambda maxes are slightly more coarsly spaced than optical densities.
psiParamsDomainList = {linspace(-2,2,10),linspace(-2,2,10),...
    linspace(-2,2,10),linspace(-2,2,10),linspace(-2,2,10),...
    linspace(-2,2,9),linspace(-2,2,9),linspace(-2,2,9)};
for i = 1:length(indDiffSds)
    psiParamsDomainList{i} = psiParamsDomainList{i}*variableParams(i)*indDiffSds(i);
    psiParamsDomainList{i} = unique(psiParamsDomainList{i});
end

% Stimulus parameter list
stimParamsDomainList = {0:0.05:1,0:0.05:1,testWls};

% Psychometric function
noiseScaleFactor = 3;
lambdaRef = 0.8;
PFSim = @(stimParams,coneParams)qpPFRMFull(stimParams,coneParams,...
    opponentVec,observer.colorDiffParams.noiseSd*noiseScaleFactor,S,...
    p1Spd,p2Spd,testSpds,testWls,lambdaRef);

% Simulated observer function, as a function of stimulus parameter
for ii = 1:length(simConeParams)
    if (simConeParams(ii) ~= 0 & variableParams(ii) == 0)
        error('Varying a simulated parameter that we are telling Quest+ is locked');
    end
end
simObserverFun = @(stimParams) qpSimulatedObserver(stimParams,PFSim,simConeParams);

% Set up a Quest object, with an option to use precomputed data if
% available
USE_PRECOMPUTE = false;
if (~USE_PRECOMPUTE)
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
questData = questDataRaw;

% Run trials
nTrials = 300;
trialPrintout = 20;
startTime = tic;
for tt = 1:nTrials
    % Get stimulus for this trial
    stim = qpQuery(questData);
    
    % Simulate outcome
    outcome = simObserverFun(stim);
    
    % Update quest data structure
    questData = qpUpdate(questData,stim,outcome);
    
    % Periodic prientout
    if (rem(tt,trialPrintout) == 0)
        fprintf('\tTrial %d of %d\n',tt,nTrials);
    end
end
elapsedTime = toc(startTime);
fprintf('Done with trial simulation, %0.3f calculation time per trial\n',(elapsedTime)/nTrials);
save('qDataSim','questData');
disp('Done Looping');

% Max posterior estimate
psiParamsIndex = qpListMaxArg(questData.posterior);
psiParamsQuest = questData.psiParamsDomain(psiParamsIndex,:);
fprintf('Simulated cone parameters: %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f, %0.3f', ...
    simConeParams(1),simConeParams(2),simConeParams(3),simConeParams(4), ...
    simConeParams(5),simConeParams(6),simConeParams(7),simConeParams(8));
fprintf('\n');
fprintf('Max posterior QUEST+ cone parameters: %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f, %0.3f', ...
    psiParamsQuest(1),psiParamsQuest(2),psiParamsQuest(3),psiParamsQuest(4), ...
    psiParamsQuest(5),psiParamsQuest(6),psiParamsQuest(7),psiParamsQuest(8));
fprintf('\n');

% Maximum likelihood fit.
% Use psiParams from QUEST+ as the starting parameter for the search, and
% impose as parameter bounds the range provided to QUEST+.
[domainVlb,domainVub] = qpGetBoundsFromDomainList(psiParamsDomainList);
psiParamsFit = qpFit(questData.trialData,questData.qpPF,psiParamsQuest,...
    questData.nOutcomes,'lowerBounds', domainVlb,'upperBounds',domainVub);
fprintf('Maximum likelihood fit parameters:   %0.3f, %0.3f, %0.3f, %0.3f. %0.3f, %0.3f, %0.3f, %0.3f, %0.3f', ...
    psiParamsFit(1),psiParamsFit(2),psiParamsFit(3),psiParamsFit(4), ...
    psiParamsFit(5),psiParamsFit(6),psiParamsFit(7),psiParamsFit(8));
fprintf('\n');

%% Plotting
% Plot of QUEST runs
figure(1);
stimCounts = qpCounts(qpData(questData.trialData),questData.nOutcomes);
xVals = [];
yVals = [];
zVals = [];
markerColors = [];
markerSizes = [];
for cc = 1:length(stimCounts)
    stim = stimCounts(cc).stim;
    for jj = 1:questData.nOutcomes
        outcomeCount = stimCounts(cc).outcomeCounts(jj);
        if outcomeCount == 0 
            continue;
        end 
%         switch (jj)
%             case 1
%                 theColor = 'r';
%             case 2
%                 theColor = 'g';
%             case 3
%                 theColor = 'b';
%             case 4
%                 theColor = 'y';
%         end
        xVals = [xVals, stim(1)];
        yVals = [yVals, stim(2)];
        zVals = [zVals, stim(3)];
        markerColors = [markerColors, jj];
        markerSizes = [markerSizes,1000*outcomeCount/max(nTrials)];
    end
end
scatter3(xVals,yVals,zVals,markerSizes,markerColors,'o','LineWidth',2,...
    'MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4);
xlim([0 1]);
ylim([0 1]);
zlim([min(testWls)-10,max(testWls)+10]);
xlabel('Lambda (Proportion Red)');
ylabel('Test Intensity');
zlabel('Peak Test Wavelength (nm)');
title('QUEST+ Trial Placement');

% Plot of cone parameters
paramsFig = figure(2);
set(paramsFig,'Color',[1 1 1],'Position',[10 10 1700 800]);
hold on;
nCols2 = 4;
nRows2 = ceil(length(psiParamsFit)/nCols2);
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
for ii = 1:length(psiParamsFit)
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
    xVals = simConeParams(ii);   % Predicted parameters
    yVals = psiParamsFit(ii); % Recovered params
    plot(xVals,yVals,'b* ','MarkerSize',7,'LineWidth',1);
    refline(1,0);
    
    % Titles and labels
    theTitle = sprintf('%s Recovered vs Simulated',cell2mat(coneParamNames(ii)));
    title(theTitle);
    xlabel('Simulated Parameters');
    ylabel('Recovered Parameters');
end

% Plots of cone fundamentals
observerRecovered = genRayleighObserver('coneVec',psiParamsFit,...
    'opponentParams',opponentVec,'S',S);
figure(3)
hold on;
h1 = plot(wls,observer.T_cones(1:2,:),'r-','LineWidth',2.5);
h2 = plot(wls,observerRecovered.T_cones(1:2,:),'b-','LineWidth',1.5);
legend([h1(1) h2(1)],'Simulated Cones', 'Recovered Cones'); 
title('L and M Cones');
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

figure(4); clf;
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