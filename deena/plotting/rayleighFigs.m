% Generates figures for Rayleigh matching manuscript

% History
%     11/5/20  dce   -Wrote it
%     11/8/20  dce   -Added opponent sphere
%     12/1/20  dce   -Fixed Pitt diagram

%% Setup
% Useful parameters
S = [380 2 201];
wls = SToWls(S);
standardObs = genRayleighObserver('S',S);
coneParamSds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];

% Initialize plotLab recipe
plotlabOBJ = plotlab();
plotlabOBJ.applyRecipe(...
    'colorOrder', [1 0 0; 0 0 1; 0 1 0], ...
    'lineMarkerSize', 12, ...
    'figureWidthInches', 8, ...
    'figureHeightInches', 5,...
    'legendLocation', 'NorthEast',...
    'axesFontSize', 16,...
    'axesLabelFontSizeMultiplier', 1,...
    'axesMinorGridAlpha', 1);

%% 4 - panel figure with cone case study

% Generate observers with the desired parameters
lODParams = [0 0 2 0 0 0 0 0].*coneParamSds;
lLambdaMaxParams = [0 0 0 0 0 2 0 0].*coneParamSds;
odObs = genRayleighObserver('coneVec',lODParams,'S',S);
lambdaMaxObs = genRayleighObserver('coneVec',lLambdaMaxParams,'S',S);

% Set up figure and subplots
coneRecoveryPlotNoiseless = figure(1);
hold on;
theAxesGrid = plotlab.axesGrid(coneRecoveryPlotNoiseless, ...
    'rowsNum', 2, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Plot results
set(gcf,'CurrentAxes',theAxesGrid{1,1});
axis(theAxesGrid{1,1}, 'square');
plot(wls,standardObs.T_cones(1,:),wls,lambdaMaxObs.T_cones(1,:));
theTitle = sprintf('L cone Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Standard', 'Shifted','FontSize',12);
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
plot(wls,standardObs.T_cones(1,:),wls,odObs.T_cones(1,:));
theTitle = sprintf('L cone Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Standard', 'Shifted','FontSize',12);
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:));
xlabel('Wavelength (nm)');
ylabel('Shifted Cone - Standard Cone');
ylim([-0.06 0.06]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:));
xlabel('Wavelength (nm)');
ylabel('Shifted Cone - Standard Cone');
ylim([-0.06 0.06]);

%% Pitt diagrams
% Produces two figures. The first is a two panel figure where the first
% panel show how match patterns change with lambda max, and the second
% shows how this is obscured by optical density shifts. The second figure
% shows how the patches can be separated by using different reference wavelengths

% Initialize figure
% Set up figure and subplots
pittPlot = figure(2);
hold on;
theAxesGrid = plotlab.axesGrid(pittPlot, ...
    'rowsNum', 1, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Fix axis limits
for i = 1:2
    theAxesGrid{1,i}.XTick = [0 0.2 0.4 0.6 0.8 1];
    theAxesGrid{1,i}.YTick = [0 0.2 0.4 0.6 0.8 1];
    theAxesGrid{1,i}.XTickLabel = {'0','0.2','0.4','0.6','0.8','1'};
    theAxesGrid{1,i}.YTickLabel = {'0','0.2','0.4','0.6','0.8','1'};
end 

% Define observers
lambdaMaxParams1 = [0 0 0 0 0 2 0 0].*coneParamSds;
lambdaMaxParams2 = [0 0 0 0 0 0 0 0].*coneParamSds;
lambdaMaxODParams1 = [0 0 -2 0 0 2 0 0].*coneParamSds;
lambdaMaxODParams2 = [0 0 2 0  0 0 0 0].*coneParamSds;
opponentParams = [40.3908 205.7353 62.9590 1.0000];
S = [380 1 401];

lambdaMaxObs1 = genRayleighObserver('S',S,'opponentParams',opponentParams,...
    'coneVec',lambdaMaxParams1);
lambdaMaxObs2 = genRayleighObserver('S',S,'opponentParams',opponentParams,...
    'coneVec',lambdaMaxParams2);
lambdaMaxODObs1 = genRayleighObserver('S',S,'opponentParams',opponentParams,...
    'coneVec',lambdaMaxODParams1);
lambdaMaxODObs2 = genRayleighObserver('S',S,'opponentParams',opponentParams,...
    'coneVec',lambdaMaxODParams2);

% Define spectra
wls = SToWls(S);
p1Wl = 560;
p2Wl = 670;
testWlInitial = 590;
testWls = [570 590 600 610 620 630];
theLegend = {'570','590','600','610','620','630'};

p1Spd = zeros(length(wls),1);
p2Spd = zeros(length(wls),1);
testSpds = zeros(length(wls),length(testWls));

p1Spd(wls==p1Wl) = 1;
p2Spd(wls==p2Wl) = 1;
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1;
end

% Plot results
noiseScalar = 4;
colors = {'Red','Blue'};

set(gcf,'CurrentAxes',theAxesGrid{1,1});
theTitle = 'Vary L $\lambda_{max}$';
hold on;
plotRayleighMatchesObserver(lambdaMaxObs1,p1Spd,p2Spd,testSpds(:,testWls==testWlInitial),...
    noiseScalar,colors{1},theTitle,'figHandle',pittPlot)
plotRayleighMatchesObserver(lambdaMaxObs2,p1Spd,p2Spd,testSpds(:,testWls==testWlInitial),...
    noiseScalar,colors{2},theTitle,'figHandle',pittPlot)
legend('$\lambda_{max}$ +2 Sd','Standard Observer','interpreter','latex')
ylabel('Reference Intensity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
theTitle = 'Vary L $\lambda_{max}$ and Optical Density';
plotRayleighMatchesObserver(lambdaMaxODObs1,p1Spd,p2Spd,testSpds(:,testWls==testWlInitial),...
    noiseScalar,colors{1},theTitle,'figHandle',pittPlot)
plotRayleighMatchesObserver(lambdaMaxODObs2,p1Spd,p2Spd,testSpds(:,testWls==testWlInitial),...
    noiseScalar,colors{2},theTitle,'figHandle',pittPlot)
legend('$\lambda_{max}$ +2 Sd, OD -2 Sd','OD +2 Sd','interpreter','latex')
ylabel('Reference Intensity');

gPittPlot = figure(3);
theTitle = 'Generalized Pitt Diagram: Vary L $\lambda_{max}$ and Optical Density';
for i = 1:length(testWls)
    plotRayleighMatchesObserver(lambdaMaxODObs1,p1Spd,p2Spd,testSpds(:,i),...
        noiseScalar,colors{1},theTitle,'figHandle',gPittPlot)
    plotRayleighMatchesObserver(lambdaMaxODObs2,p1Spd,p2Spd,testSpds(:,i),...
        noiseScalar,colors{2},theTitle,'figHandle',gPittPlot)
end
legend('$\lambda_{max}$ +2 Sd, OD -2 Sd','OD +2 Sd','interpreter','latex')
ylabel('Reference Intensity');

%% Opponent contrast sphere figure
noiseSD = 1;
[~,~,opponentSphere] = opponentAxesToLab(noiseSD);
opponentSphereFig = figure(4);
plot3(opponentSphere(1,:),opponentSphere(2,:),opponentSphere(3,:),'o-');
xlabel('x');
ylabel('y');
zlabel('z');
title('Optimized Opponent Sphere, Radius = 1');

%% Peak spectra
lightFile = 'OLRayleighMatch3201SpectralSettings_670_560_600_1_0.02_0.5.mat';
lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops',lightFile);
lightData = load(lightFileName);
S = [380 2 201];
SNew = [380 1 401];
wls = SToWls(SNew);

peakSpectralPlot = figure(5);
hold on;
theAxesGrid = plotlab.axesGrid(peakSpectralPlot, ...
    'rowsNum', 1, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

set(gcf,'CurrentAxes',theAxesGrid{1,1});
hold on;
p1Spd = SplineSpd(S,lightData.primary1IncrSpd,SNew,0);
p2Spd = SplineSpd(S,lightData.primary2IncrSpd,SNew,0);
plot(wls,p1Spd.*lightData.p1ScaleFactor,'r');
plot(wls,p2Spd*lightData.p2ScaleFactor,'g');
legend('Primary 1','Primary 2','Location', 'NorthWest');
title('Primary Spds');
xlabel('Wavelength (nm)');
ylabel('Radiance (W/[sr-m^2-nm])');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
testSpd = SplineSpd(S,lightData.testIncrSpd,SNew,0);
plot(wls,testSpd*lightData.testScaleFactor,'y');
legend('Reference','Location', 'NorthWest');
title('Reference Spd (600 nm)');
xlabel('Wavelength (nm)');
ylabel('Radiance (W/[sr-m^2-nm])');

%% New recipe
plotlabOBJ = plotlab();
plotlabOBJ.applyRecipe(...
    'colorOrder', [1 0 0; 0 0 1; 0 1 0], ...
    'lineMarkerSize', 12, ...
    'figureWidthInches', 8, ...
    'figureHeightInches', 8,...
    'legendLocation', 'NorthEast',...
    'axesFontSize', 12,...
    'axesLabelFontSizeMultiplier', 1,...
    'axesMinorGridAlpha', 1);

%% Rayleigh Simulation Plot - Noiseless Cone Recovery Figure
% Generate observers with the desired parameters
lODParams = [0 0 2 0 0 0 0 0].*coneParamSds;
lLambdaMaxParams = [0 0 0 0 0 2 0 0].*coneParamSds;
S = [380 2 201];
opponentParams = [40.3908 205.7353 62.9590 1.0000];
odObs = genRayleighObserver('coneVec',lODParams,'S',S,'opponentParams',...
    opponentParams);
lambdaMaxObs = genRayleighObserver('coneVec',lLambdaMaxParams,'S',S,...
    'opponentParams',opponentParams);

% Recover parameters

resFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','coneShiftNoiseless','coneShiftNoiseless_paramsSearchData');
coneShiftData = load(resFile);
recoveredParams = coneShiftData.recoveredParams;

odObsRec = genRayleighObserver('coneVec',recoveredParams(1,:),'S',S,...
    'opponentParams',opponentParams);
lambdaMaxObsRec = genRayleighObserver('coneVec',recoveredParams(2,:),'S',S,...
    'opponentParams',opponentParams);

% Set up figure and subplots
coneRecoveryPlotNoisy = figure(7);
hold on;
theAxesGrid = plotlab.axesGrid(coneRecoveryPlotNoisy, ...
    'rowsNum', 2, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Plot results
set(gcf,'CurrentAxes',theAxesGrid{1,1});
axis(theAxesGrid{1,1}, 'square');
hold on;
plot(wls,standardObs.T_cones(1,:),'g','LineWidth',4);
plot(wls,lambdaMaxObs.T_cones(1,:),'b','LineWidth',2.5);
plot(wls,lambdaMaxObsRec.T_cones(1,:),'r','LineWidth',1);
theTitle = sprintf('L Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
hold on;
plot(wls,standardObs.T_cones(1,:),'g','LineWidth',4);
plot(wls,odObs.T_cones(1,:),'b','LineWidth',2.5);
plot(wls,odObsRec.T_cones(1,:),'r','LineWidth',1);
theTitle = sprintf('L Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
hold on;
plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3)
plot(wls,lambdaMaxObs.T_cones(1,:)-lambdaMaxObsRec.T_cones(1,:),'r','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Sensitivity Difference');
lgd = legend('Simulated - Standard', 'Simulated - Recovered','Location','south');
ylim([-0.085 0.085]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
hold on;
plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
plot(wls,odObs.T_cones(1,:)-odObsRec.T_cones(1,:),'r','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Sensitivity Difference');
lgd = legend('Simulated - Standard', 'Simulated - Recovered','Location','south');
ylim([-0.085 0.085]);

sgtitle('Cone Recovery - Noiseless Case','FontSize',18);

%% Rayleigh Simulation Plot - Noisy Cone Recovery Figure
lODParams = [0 0 2 0 0 0 0 0].*coneParamSds;
lLambdaMaxParams = [0 0 0 0 0 2 0 0].*coneParamSds;
S = [380 2 201];
opponentParams = [40.3908 205.7353 62.9590 1.0000];
odObs = genRayleighObserver('coneVec',lODParams,'S',S,'opponentParams',...
    opponentParams);
lambdaMaxObs = genRayleighObserver('coneVec',lLambdaMaxParams,'S',S,...
    'opponentParams',opponentParams);

% Recover parameters

resFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','coneShiftNoisy','coneShiftNoisy_paramsSearchData');
coneShiftData = load(resFile);
recoveredParams = coneShiftData.recoveredParams;

odObsRec = genRayleighObserver('coneVec',recoveredParams(1,:),'S',S,...
    'opponentParams',opponentParams);
lambdaMaxObsRec = genRayleighObserver('coneVec',recoveredParams(2,:),'S',S,...
    'opponentParams',opponentParams);

% Set up figure and subplots
coneRecoveryPlotNoiseless = figure(6);
hold on;
theAxesGrid = plotlab.axesGrid(coneRecoveryPlotNoiseless, ...
    'rowsNum', 2, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Plot results
set(gcf,'CurrentAxes',theAxesGrid{1,1});
axis(theAxesGrid{1,1}, 'square');
hold on;
plot(wls,standardObs.T_cones(1,:),'g','LineWidth',4);
plot(wls,lambdaMaxObs.T_cones(1,:),'b','LineWidth',2.5);
plot(wls,lambdaMaxObsRec.T_cones(1,:),'r','LineWidth',1);
theTitle = sprintf('L Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)-0.01])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
hold on;
plot(wls,standardObs.T_cones(1,:),'g','LineWidth',4);
plot(wls,odObs.T_cones(1,:),'b','LineWidth',2.5);
plot(wls,odObsRec.T_cones(1,:),'r','LineWidth',1);
theTitle = sprintf('L Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)-0.01])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
hold on;
plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3)
plot(wls,lambdaMaxObs.T_cones(1,:)-lambdaMaxObsRec.T_cones(1,:),'r','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Sensitivity Difference');
lgd = legend('Simulated - Standard', 'Simulated - Recovered','Location','south');
ylim([-0.085 0.085]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
hold on;
plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
plot(wls,odObs.T_cones(1,:)-odObsRec.T_cones(1,:),'r','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Sensitivity Difference');
lgd = legend('Simulated - Standard', 'Simulated - Recovered','Location','south');
ylim([-0.085 0.085]);

sgtitle('Cone Recovery: Noise Sd = 4','FontSize',18);

%% Parameter recovery plots
% Noiseless case
paramRecoveryPlotNoiseless = figure(8);
hold on;
theAxesGrid = plotlab.axesGrid(paramRecoveryPlotNoiseless, ...
    'rowsNum', 2, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

noiselessParamFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper2_0_FC','varyNoiseForPaper2_0_FC_paramsSearchData.mat');
paramData_0 = load(noiselessParamFile);
coneParamNames = {'L Optical Density','M Optical Density',...
    'L Lambda Max','M Lambda Max'};
nCols = 2;
paramInds = [3 4 6 7]; % Variable params, to plot

for ii = 1:4
    % Make a subplot in the correct position
    row = ceil(ii/nCols);
    col = mod(ii,nCols);
    if col == 0
        col = nCols;
    end
    set(gcf,'CurrentAxes',theAxesGrid{row,col});
    axis(theAxesGrid{row,col}, 'square');
    hold on;
    
    % Define axis limits
    if ii < 3  % Density shifts, in %
        limits = [-40 40];
    else        % Wavelength shifts, in nm
        limits = [-5 5];
    end
    xlim(limits);
    ylim(limits);
    
    % Set titles and labels
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
    title(cell2mat(coneParamNames(ii)));
    
    % Plot data 
    xVals = paramData_0.sampledParams(:,paramInds(ii));   % Predicted parameters
    yVals = paramData_0.recoveredParams(:,paramInds(ii)); % Recovered params
    l2 = refline(1,0);
    l1 = plot(xVals,yVals,'. ','MarkerEdgeColor','Blue');
end
sgtitle('Recovered vs Simulated Parameters: Noiseless Case','FontSize',18);

%% New recipe
% Initialize plotLab recipe
plotlabOBJ = plotlab();
plotlabOBJ.applyRecipe(...
    'colorOrder', [1 0 0; 0 0 1; 0 1 0], ...
    'lineMarkerSize', 12, ...
    'figureWidthInches', 8, ...
    'figureHeightInches', 5,...
    'legendLocation', 'NorthEast',...
    'axesFontSize', 12,...
    'axesLabelFontSizeMultiplier', 1,...
    'axesMinorGridAlpha', 1);

%% Noisy case 
paramData_0 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper2_0_FC','varyNoiseForPaper2_0_FC_paramsSearchData.mat'));
paramData_1 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper2_1_FC','varyNoiseForPaper2_1_FC_paramsSearchData.mat'));
paramData_2 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper2_2_FC','varyNoiseForPaper2_2_FC_paramsSearchData.mat'));
paramData_3 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper2_3_FC','varyNoiseForPaper2_3_FC_paramsSearchData.mat'));
paramData_4 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper2_4_FC','varyNoiseForPaper2_4_FC_paramsSearchData.mat'));
noisyData = {paramData_0 paramData_1 paramData_2 paramData_3 paramData_4};

odRecoveryPlot = figure(9);
hold on;
theAxesGridOD = plotlab.axesGrid(odRecoveryPlot, ...
    'rowsNum', 2, 'colsNum', 3, ...
    'spacing', 'compact', ...
    'padding', 'compact', ...
    'method', 'tile');

lambdaMaxRecoveryPlot = figure(10);
hold on;
theAxesGridLambda = plotlab.axesGrid(lambdaMaxRecoveryPlot, ...
    'rowsNum', 2, 'colsNum', 3, ...
    'spacing', 'compact', ...
    'padding', 'compact', ...
    'method', 'tile');

nRow = 2;
nCol = 3;
for i = 1:5
    row = ceil(i/nCol); 
    col = mod(i,nCol);
    if col == 0
        col = nCol;
    end
    theTitle = sprintf('Noise Sd = %s',num2str(i)-1);
    
    figure(9)
    set(gcf,'CurrentAxes',theAxesGridOD{row,col});
    axis(theAxesGridOD{row,col}, 'square');
    hold on;
    limits = [-40 40];
    xlim(limits);
    ylim(limits);
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
    title(theTitle);
    refline(1,0);
    l1 = plot(noisyData{i}.sampledParams(:,3),noisyData{i}.recoveredParams(:,3),...
        '. ','MarkerEdgeColor','Yellow');
    l2 = plot(noisyData{i}.sampledParams(:,4),noisyData{i}.recoveredParams(:,4),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
    
    figure(10)
    set(gcf,'CurrentAxes',theAxesGridLambda{row,col});
    axis(theAxesGridLambda{row,col},'square');
    hold on;
    limits = [-5 5];
    xlim(limits);
    ylim(limits);
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
    title(theTitle);
    refline(1,0);
    l1 = plot(noisyData{i}.sampledParams(:,6),noisyData{i}.recoveredParams(:,6),...
        '. ','MarkerEdgeColor','Yellow');
    l2 = plot(noisyData{i}.sampledParams(:,7),noisyData{i}.recoveredParams(:,7),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
end 
figure(9);
sgtitle('Optical Density Recovered vs Simulated Parameters: Noisy','FontSize',18);

figure(10);
sgtitle('Lambda Max Recovered vs Simulated Parameters: Noisy','FontSize',18);

%% Noisy case - vary number of matches per test wavelength 

% Initialize plotLab recipe
plotlabOBJ = plotlab();
plotlabOBJ.applyRecipe(...
    'colorOrder', [1 0 0; 0 0 1; 0 1 0], ...
    'lineMarkerSize', 12, ...
    'figureWidthInches', 9, ...
    'figureHeightInches', 8,...
    'legendLocation', 'NorthEast',...
    'axesFontSize', 12,...
    'axesLabelFontSizeMultiplier', 1,...
    'axesMinorGridAlpha', 1);

paramData_1Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper7_1_FC','varyNMatchesForPaper7_1_FC_paramsSearchData.mat'));
paramData_2Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper7_2_FC','varyNMatchesForPaper7_2_FC_paramsSearchData.mat'));
paramData_3Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper7_3_FC','varyNMatchesForPaper7_3_FC_paramsSearchData.mat'));
paramData_4Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper7_4_FC','varyNMatchesForPaper7_4_FC_paramsSearchData.mat'));
nMatchData = {paramData_1Match paramData_2Match paramData_3Match paramData_4Match};

odRecoveryPlot = figure(11);
hold on;
theAxesGridOD = plotlab.axesGrid(odRecoveryPlot, ...
    'rowsNum', 2, 'colsNum', 2, ...
    'spacing', 'compact', ...
    'padding', 'compact', ...
    'method', 'tile');

lambdaMaxRecoveryPlot = figure(12);
hold on;
theAxesGridLambda = plotlab.axesGrid(lambdaMaxRecoveryPlot, ...
    'rowsNum', 2, 'colsNum', 2, ...
    'spacing', 'compact', ...
    'padding', 'compact', ...
    'method', 'tile');

nRow = 2;
nCol = 2;
for i = 1:4
    row = ceil(i/nCol); 
    col = mod(i,nCol);
    if col == 0
        col = nCol;
    end
    theTitle = sprintf('%s Matches per Reference Wavelength',num2str(i));
    
    figure(11)
    set(gcf,'CurrentAxes',theAxesGridOD{row,col});
    axis(theAxesGridOD{row,col}, 'square');
    hold on;
    limits = [-40 40];
    xlim(limits);
    ylim(limits);
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
    title(theTitle);
    refline(1,0);
    l1 = plot(nMatchData{i}.sampledParams(:,3),nMatchData{i}.recoveredParams(:,3),...
        '. ','MarkerEdgeColor','Yellow');
    l2 = plot(nMatchData{i}.sampledParams(:,4),nMatchData{i}.recoveredParams(:,4),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
    
    figure(12)
    set(gcf,'CurrentAxes',theAxesGridLambda{row,col});
    axis(theAxesGridLambda{row,col},'square');
    hold on;
    limits = [-5 5];
    xlim(limits);
    ylim(limits);
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
    title(theTitle);
    refline(1,0);
    l1 = plot(nMatchData{i}.sampledParams(:,6),nMatchData{i}.recoveredParams(:,6),...
        '. ','MarkerEdgeColor','Yellow');
    l2 = plot(nMatchData{i}.sampledParams(:,7),nMatchData{i}.recoveredParams(:,7),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
end 
figure(11);
sgtitle('Optical Density Recovered vs Simulated: Vary Number of Matches','FontSize',18);

figure(12);
sgtitle('Lambda Max Recovered vs Simulated: Vary Number of Matches','FontSize',18);