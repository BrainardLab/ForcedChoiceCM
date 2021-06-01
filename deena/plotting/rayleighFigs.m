% Generates figures for Rayleigh matching manuscript

% History
%     11/5/20  dce   -Wrote it
%     11/8/20  dce   -Added opponent sphere
%     12/1/20  dce   -Fixed Pitt diagram
%     4/30/21  dce   -Changed plot colors
%     5/08/21  dce   -New figures
%     5/10/21  dce   -Updated which simulation data files we use

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
    'figureHeightInches', 8,...
    'legendLocation', 'NorthEast',...
    'axesFontSize', 16,...
    'axesLabelFontSizeMultiplier', 1,...
    'axesMinorGridAlpha', 1);
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
hold on;
plot(wls,standardObs.T_cones(1,:),'r','LineWidth',5)
plot(wls,lambdaMaxObs.T_cones(1,:),'g','LineWidth',2);
theTitle = sprintf('L cone Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated','FontSize',12);
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
hold on;
plot(wls,standardObs.T_cones(1,:),'r','LineWidth',5)
plot(wls,odObs.T_cones(1,:),'g','LineWidth',2);
theTitle = sprintf('L cone Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated','FontSize',12);
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
hold on;
plot(wls,zeros(1,length(wls)),'r--','LineWidth',6);
a = plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
legend(a, 'Simulated - Standard');
ylim([-0.085 0.085]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
hold on;
plot(wls,zeros(1,length(wls)),'r--','LineWidth',6);
a = plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
legend(a, 'Simulated - Standard');
ylim([-0.085 0.085]);

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
coneRecoveryPlotNoiseless = figure(7);
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
plot(wls,lambdaMaxObs.T_cones(1,:),'g','LineWidth',5);
plot(wls,lambdaMaxObsRec.T_cones(1,:),'b','LineWidth',2);
theTitle = sprintf('L Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
hold on;
plot(wls,odObs.T_cones(1,:),'g','LineWidth',5);
plot(wls,odObsRec.T_cones(1,:),'b','LineWidth',2);
theTitle = sprintf('L Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
hold on;
plot(wls,zeros(1,length(wls)),'r--','LineWidth',6);
a = plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
b = plot(wls,lambdaMaxObs.T_cones(1,:)-lambdaMaxObsRec.T_cones(1,:),'b','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
lgd = legend([a b],'Simulated - Standard', 'Simulated - Recovered','Location','south');
ylim([-0.085 0.085]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
hold on;
plot(wls,zeros(1,length(wls)),'r--','LineWidth',6);
a = plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
b = plot(wls,odObs.T_cones(1,:)-odObsRec.T_cones(1,:),'b','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
lgd = legend([a b],'Simulated - Standard', 'Simulated - Recovered','Location','south');
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

odObsRecNoisy = genRayleighObserver('coneVec',recoveredParams(1,:),'S',S,...
    'opponentParams',opponentParams);
lambdaMaxObsRecNoisy = genRayleighObserver('coneVec',recoveredParams(2,:),'S',S,...
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
plot(wls,lambdaMaxObs.T_cones(1,:),'g','LineWidth',5);
plot(wls,lambdaMaxObsRecNoisy.T_cones(1,:),'b','LineWidth',2);
theTitle = sprintf('L Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)-0.01])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
hold on;
plot(wls,odObs.T_cones(1,:),'g','LineWidth',5);
plot(wls,odObsRecNoisy.T_cones(1,:),'b','LineWidth',2);
theTitle = sprintf('L Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Simulated', 'Recovered');
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)-0.01])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
hold on;
plot(wls,zeros(1,length(wls)),'r--','LineWidth',6);
a = plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
b = plot(wls,lambdaMaxObs.T_cones(1,:)-lambdaMaxObsRecNoisy.T_cones(1,:),'b','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Sensitivity Difference');
lgd = legend([a b],'Simulated - Standard', ...
    'Simulated - Recovered','Location','south');
ylim([-0.085 0.085]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
hold on;
plot(wls,zeros(1,length(wls)),'r--','LineWidth',6);
a = plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:),'g','LineWidth',3);
b = plot(wls,odObs.T_cones(1,:)-odObsRecNoisy.T_cones(1,:),'b','LineWidth',3);
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Sensitivity Difference');
lgd = legend([a b],'Simulated - Standard',...
    'Simulated - Recovered','Location','south');
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
    'paramsSearch','varyNoiseForPaper10_0_FC','varyNoiseForPaper10_0_FC_paramsSearchData.mat');
paramData_0 = load(noiselessParamFile);
coneParamNames = {'L Lambda Max','M Lambda Max','L Optical Density','M Optical Density'};
nCols = 2;
paramInds = [6 7 3 4]; % Variable params, to plot

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
    if ii >= 3  % Density shifts, in %
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
    l2 = plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
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
    'paramsSearch','varyNoiseForPaper10_0_FC','varyNoiseForPaper10_0_FC_paramsSearchData.mat'));
paramData_1 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_1_FC','varyNoiseForPaper10_1_FC_paramsSearchData.mat'));
paramData_2 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_2_FC','varyNoiseForPaper10_2_FC_paramsSearchData.mat'));
paramData_3 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_3_FC','varyNoiseForPaper10_3_FC_paramsSearchData.mat'));
paramData_4 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_4_FC','varyNoiseForPaper10_4_FC_paramsSearchData.mat'));
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
    plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
    l1 = plot(noisyData{i}.sampledParams(:,3),noisyData{i}.recoveredParams(:,3),...
        '. ','MarkerEdgeColor','Cyan');
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
    plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
    l1 = plot(noisyData{i}.sampledParams(:,6),noisyData{i}.recoveredParams(:,6),...
        '. ','MarkerEdgeColor','Cyan');
    l2 = plot(noisyData{i}.sampledParams(:,7),noisyData{i}.recoveredParams(:,7),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
end 
figure(9);
sgtitle('Optical Density Recovered vs Simulated Parameters: Noisy','FontSize',18);

figure(10);
sgtitle('Lambda Max Recovered vs Simulated Parameters: Noisy','FontSize',18);



%% Noisy case- cone excitations fitting method
% Note that this fit is based on cone excitations and is unaveraged
paramData_0 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_0_FC','varyNoiseForPaper10_0_FC_paramsSearchData_ExcitationFits3.mat'));
paramData_1 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_1_FC','varyNoiseForPaper10_1_FC_paramsSearchData_ExcitationFits3.mat'));
paramData_2 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_2_FC','varyNoiseForPaper10_2_FC_paramsSearchData_ExcitationFits3.mat'));
paramData_3 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_3_FC','varyNoiseForPaper10_3_FC_paramsSearchData_ExcitationFits3.mat'));
paramData_4 = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNoiseForPaper10_4_FC','varyNoiseForPaper10_4_FC_paramsSearchData_ExcitationFits3.mat'));
noisyData = {paramData_0 paramData_1 paramData_2 paramData_3 paramData_4};

ODConeRecoveryPlot = figure(11);
hold on;
theAxesGridOD = plotlab.axesGrid(ODConeRecoveryPlot, ...
    'rowsNum', 2, 'colsNum', 3, ...
    'spacing', 'compact', ...
    'padding', 'compact', ...
    'method', 'tile');

lambdaMaxConeRecoveryPlot = figure(12);
hold on;
theAxesGridLambda = plotlab.axesGrid(lambdaMaxConeRecoveryPlot, ...
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
    plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
    l1 = plot(noisyData{i}.sampledParams(:,3),noisyData{i}.recoveredParams(:,3),...
        '. ','MarkerEdgeColor','Cyan');
    l2 = plot(noisyData{i}.sampledParams(:,4),noisyData{i}.recoveredParams(:,4),...
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
    plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
    l1 = plot(noisyData{i}.sampledParams(:,6),noisyData{i}.recoveredParams(:,6),...
        '. ','MarkerEdgeColor','Cyan');
    l2 = plot(noisyData{i}.sampledParams(:,7),noisyData{i}.recoveredParams(:,7),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
end 
figure(11);
sgtitle('Optical Density Recovered vs Simulated Parameters: Cone Fit','FontSize',18);

figure(12);
sgtitle('Lambda Max Recovered vs Simulated Parameters: Cone Fit','FontSize',18);

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
    'paramsSearch','varyNMatchesForPaper11_1_FC','varyNMatchesForPaper11_1_FC_paramsSearchData.mat'));
paramData_2Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper11_2_FC','varyNMatchesForPaper11_2_FC_paramsSearchData.mat'));
paramData_3Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper11_3_FC','varyNMatchesForPaper11_3_FC_paramsSearchData.mat'));
paramData_4Match = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearch','varyNMatchesForPaper11_4_FC','varyNMatchesForPaper11_4_FC_paramsSearchData.mat'));
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
    plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
    l1 = plot(nMatchData{i}.sampledParams(:,3),nMatchData{i}.recoveredParams(:,3),...
        '. ','MarkerEdgeColor','Cyan');
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
    plot(limits(1):limits(2),limits(1):limits(2),'r--','LineWidth',2);
    l1 = plot(nMatchData{i}.sampledParams(:,6),nMatchData{i}.recoveredParams(:,6),...
        '. ','MarkerEdgeColor','Cyan');
    l2 = plot(nMatchData{i}.sampledParams(:,7),nMatchData{i}.recoveredParams(:,7),...
        '. ','MarkerEdgeColor','Blue');
    legend([l1 l2],'L','M','Location','southeast');
end 
figure(11);
sgtitle('Optical Density Recovery: Vary Number of Matches (Noise SD = 3)','FontSize',18,'Interpreter','none');

figure(12);
sgtitle('Lambda Max Recovery: Vary Number of Matches (Noise SD = 3)','FontSize',18,'Interpreter','none');


%% Human data - cone excitation ratio plots (Standard)
% Load data 
data3044 = load(fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'MELA_3044','MELA_3044_11_26_analysis_unavgNew.mat'));
data3045 = load(fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'MELA_3045','MELA_3045_11_26_analysis_unavgNew.mat'));

a = figure(13);
hold on;
b = figure(14);
hold on;

humanData = {data3044, data3045};
figHandles = {13 14};
subjIDs = {'MELA_3044', 'MELA_3045'};
plotColors = {'Red','Black','Blue','Green','Cyan','Magenta'};

% Loop through the two human observers, plot each of their data
for i = 1:2
    theData = humanData{i};
    matchWls = unique(theData.lightCombos,'rows');
    legendHandles = [];
    legendEntries = {};
    
    % Plot setup 
    figure(figHandles{i});
    xlim([0.05 0.35]);
    ylim([0.015 0.035]);
    xlabel('(L - M)/(L+M)');
    ylabel('L + M');
    title([subjIDs{i} ' Cone Response Difference - Standard Cones'],'interpreter','none');
    
    for j = 1:size(matchWls,1) % Loop through each wavelength
        % Extract relevant match spds
        sessionInds = all(theData.lightCombos==matchWls(j,:),2);      
        measPrimarySpdsTrial = theData.measPrimarySpds(:,sessionInds);
        measRefSpdsTrial = theData.measRefSpds(:,sessionInds);
        
        % Compute excitations
        primaryResStd = theData.stdObs.T_cones*measPrimarySpdsTrial;
        refResStd = theData.stdObs.T_cones*measRefSpdsTrial;
        primaryLMinusMStd = primaryResStd(1,:)-primaryResStd(2,:);
        primaryLPlusMStd = primaryResStd(1,:)+primaryResStd(2,:);
        refLMinusMStd = refResStd(1,:)-refResStd(2,:);
        refLPlusMStd = refResStd(1,:)+refResStd(2,:);
        
        % Add data to plot 
        line1  = plot(primaryLMinusMStd./primaryLPlusMStd,...
            primaryLPlusMStd,'* ','MarkerEdgeColor',...
            plotColors{j},'MarkerSize',7);
        plot(refLMinusMStd./refLPlusMStd,refLPlusMStd,'o ','MarkerEdgeColor',...
            plotColors{j},'MarkerFaceColor',plotColors{j},'MarkerSize',7);
        plot(refLMinusMStd(3:4)./refLPlusMStd(3:4),refLPlusMStd(3:4),...
            'yo','MarkerFaceColor','Yellow','MarkerEdgeColor',plotColors{j},...
            'MarkerSize',7);
        legendHandles = [legendHandles,line1];
        legendEntries{end+1} = num2str(matchWls(j,3));
        
        % Add lines connecting each pair of primary/reference points
        for k = 1:length(primaryResStd(1,:))
            plot([primaryLMinusMStd(k)/primaryLPlusMStd(k)...
                refLMinusMStd(k)/refLPlusMStd(k)],...
                [primaryLPlusMStd(k) refLPlusMStd(k)],'k-');
        end
    end
    % Add legend 
    legend(legendHandles,legendEntries);  
end

%% Human data - cone excitation ratio plots (fits)
% Load data 
data3044_OC = load(fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'MELA_3044','MELA_3044_11_26_analysis_unavgNew.mat'));
data3044_Cone = load(fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'MELA_3044','MELA_3044_11_26_analysis_coneErrUnavg.mat'));
data3045_OC = load(fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'MELA_3045','MELA_3045_11_26_analysis_unavgNew.mat'));
data3045_Cone = load(fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'MELA_3045','MELA_3045_11_26_analysis_coneErrUnavg.mat'));

% Set up figures
ocFit3044 = figure(15);
hold on;
coneFit3044 = figure(16);
hold on;
ocFit3045 = figure(17);
hold on;
coneFit3045 = figure(18);
hold on;

humanData = {data3044_OC, data3044_Cone, data3045_OC, data3045_Cone};
figIDs = [15 16 17 18];
subjIDs = {'MELA_3044 Opponent Contrast Fit','MELA_3044 Cone Excitations Fit',...
    'MELA_3045 Opponent Contrast Fit',...
    'MELA_3045 Cone Excitations Fit'};
plotColors = {'Red','Black','Blue','Green','Cyan','Magenta'};

% Loop through the observers, plot their data
for i = 1:4
    theData = humanData{i};
    matchWls = unique(theData.lightCombos,'rows');
    legendHandles = [];
    legendEntries = {};
    
    % Plot setup 
    figure(figIDs(i));
    xlim([0.05 0.35]);
    ylim([0.015 0.035]);
    xlabel('(L - M)/(L+M)');
    ylabel('L + M');
    title([subjIDs{i} '-Cone Response Difference'],'interpreter','none');
    
    for j = 1:size(matchWls,1) % Loop through each wavelength
        % Extract relevant match spds
        sessionInds = all(theData.lightCombos==matchWls(j,:),2);      
        measPrimarySpdsTrial = theData.measPrimarySpds(:,sessionInds);
        measRefSpdsTrial = theData.measRefSpds(:,sessionInds);
        
        % Compute excitations
        primaryRes = theData.estObs.T_cones*measPrimarySpdsTrial;
        refRes = theData.estObs.T_cones*measRefSpdsTrial;
        primaryLMinusM = primaryRes(1,:)-primaryRes(2,:);
        primaryLPlusM = primaryRes(1,:)+primaryRes(2,:);
        refLMinusM = refRes(1,:)-refRes(2,:);
        refLPlusM = refRes(1,:)+refRes(2,:);
        
        % Add data to plot 
        line1  = plot(primaryLMinusM./primaryLPlusM,...
            primaryLPlusM,'* ','MarkerEdgeColor',...
            plotColors{j},'MarkerSize',7);
        plot(refLMinusM./refLPlusM,refLPlusM,'o ','MarkerEdgeColor',...
            plotColors{j},'MarkerFaceColor',plotColors{j},'MarkerSize',7);
        plot(refLMinusM(3:4)./refLPlusM(3:4),refLPlusM(3:4),...
            'yo','MarkerFaceColor','Yellow','MarkerEdgeColor',plotColors{j},...
            'MarkerSize',7);
        legendHandles = [legendHandles,line1];
        legendEntries{end+1} = num2str(matchWls(j,3));
        
        % Add lines connecting each pair of primary/reference points
        for k = 1:length(primaryRes(1,:))
            plot([primaryLMinusM(k)/primaryLPlusM(k)...
                refLMinusM(k)/refLPlusM(k)],...
                [primaryLPlusM(k) refLPlusM(k)],'k-');
        end
    end
    % Add legend 
    legend(legendHandles,legendEntries);  
end

%% Human Data - Cone Excitation Curve Plots
% Define directory 
baseDir = '/Users/deena/Dropbox (Aguirre-Brainard Lab)';
setpref('ForcedChoiceCM','rayleighDataDir',fullfile(baseDir,'MELA_datadev',...
    'Experiments','ForcedChoiceCM','OLRayleighMatch'));
% Set up our observers and data arrays
stdObs3044 = data3044_OC.stdObs;
coneObs3044 = data3044_Cone.estObs;
ocObs3044 = data3044_OC.estObs;
stdObs3045 = data3045_OC.stdObs;
coneObs3045 = data3045_Cone.estObs;
ocObs3045 = data3045_OC.estObs;
observers = {stdObs3044 coneObs3044 ocObs3044 stdObs3045 coneObs3045 ocObs3045};

% Each follows the columns primary L+M, primary L-M, ref L+M, ref L-M
stdObs3044Matches = [];
coneObs3044Matches = [];
ocObs3044Matches = [];
stdObs3045Matches = [];
coneObs3045Matches = [];
ocObs3045Matches = [];
coneResArrs = {stdObs3044Matches coneObs3044Matches ocObs3044Matches...
    stdObs3045Matches coneObs3045Matches ocObs3045Matches};

% Loop through possible reference wavelengths 
wls = SToWls([380 2 201]);
wls = wls(wls>570 & wls<640);
for i = 1:length(wls)
    % Find light scalars
    if wls(i) >= 630
        p2Scale = 0.004;
        testScale = 0.25;
    elseif wls(i) < 620
        p2Scale = 0.02;
        testScale = 0.1;
    else
        p2Scale = 0.004;
        testScale = 0.1;
    end 
    % Identify and load light file 
    lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
        201,670,560,wls(i),1,p2Scale,testScale);
    lightData = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops',lightFile));
    
    % Search for best match for each observer, and compute cone responses
    for j = 1:6
        [tSpd,pSpd] = searchPredictedRayleighMatch(lightData.testSpdsPredicted,...
            lightData.primarySpdsPredicted,observers{j});        
        primaryRes = theData.estObs.T_cones*pSpd;
        refRes = theData.estObs.T_cones*tSpd;
        primaryLMinusM = primaryRes(1)-primaryRes(2);
        primaryLPlusM = primaryRes(1)+primaryRes(2);
        refLMinusM = refRes(1)-refRes(2);
        refLPlusM = refRes(1)+refRes(2);
        coneResArrs{j} = [coneResArrs{j}; primaryLPlusM, primaryLMinusM,refLPlusM,refLMinusM];  
    end 
end 

% Plot results
% MELA_3044 
coneExcitationsPlot = figure(19);
theAxesGrid = plotlab.axesGrid(coneExcitationsPlot, ...
    'rowsNum', 2, 'colsNum', 1, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'topMargin',  0.03,...
    'method', 'tile');

set(gcf,'CurrentAxes',theAxesGrid{1,1});
hold on;
xlim([0.05 0.35]);
ylim([0.015 0.04]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title('MELA_3044','Interpreter','none');
a = plot(coneResArrs{1}(:,2)./coneResArrs{1}(:,1),coneResArrs{1}(:,1),'r-',...
    coneResArrs{1}(:,4)./coneResArrs{1}(:,3),coneResArrs{1}(:,3),'r--',...
    'LineWidth',6);
b = plot(coneResArrs{2}(:,2)./coneResArrs{2}(:,1),coneResArrs{2}(:,1),'g-',...
    coneResArrs{2}(:,4)./coneResArrs{2}(:,3),coneResArrs{2}(:,3),'g--',...
'LineWidth',4);
c = plot(coneResArrs{3}(:,2)./coneResArrs{3}(:,1),coneResArrs{3}(:,1),'b-',...
    coneResArrs{3}(:,4)./coneResArrs{3}(:,3),coneResArrs{3}(:,3),'b--',...
    'LineWidth',2.5);
legend([a(1) b(1) c(1)],'Standard', 'Fit (Cones)', 'Fit (Opponent Contrast');

% MELA_3045
set(gcf,'CurrentAxes',theAxesGrid{2,1});
hold on;
xlim([0.05 0.35]);
ylim([0.015 0.04]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title('MELA_3045','Interpreter','none');
a = plot(coneResArrs{4}(:,2)./coneResArrs{4}(:,1),coneResArrs{4}(:,1),'r-',...
    coneResArrs{4}(:,4)./coneResArrs{4}(:,3),coneResArrs{4}(:,3),'r--',...
    'LineWidth',6);
b = plot(coneResArrs{5}(:,2)./coneResArrs{5}(:,1),coneResArrs{5}(:,1),'g-',...
    coneResArrs{5}(:,4)./coneResArrs{5}(:,3),coneResArrs{5}(:,3),'g--',...
    'LineWidth',4);
c = plot(coneResArrs{6}(:,2)./coneResArrs{6}(:,1),coneResArrs{6}(:,1),'b-',...
    coneResArrs{6}(:,4)./coneResArrs{6}(:,3),coneResArrs{6}(:,3),'b--',...
    'LineWidth',2.5);
legend([a(1) b(1) c(1)],'Standard', 'Fit (Cones)', 'Fit (Opponent Contrast');

sgtitle('Cone Response Difference');

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

%% Rayleigh Simulation Plot - Cone L and M Recovery Figure (Noiseless)
% Generate observers with the desired parameters
lODParams = [0 0 2 0 0 0 0 0].*coneParamSds;
lLambdaMaxParams = [0 0 0 0 0 2 0 0].*coneParamSds;
S = [380 2 201];
opponentParams = [40.3908 205.7353 62.9590 1.0000];
odObs = genRayleighObserver('coneVec',lODParams,'S',S,'opponentParams',...
    opponentParams);
lambdaMaxObs = genRayleighObserver('coneVec',lLambdaMaxParams,'S',S,...
    'opponentParams',opponentParams);
stdObs = genRayleighObserver('coneVec',zeros(1,8),'S',S,...
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

% Find spds
lambdaMaxObsPrimarySpds = [];
lambdaMaxObsRefSpds = [];
odObsPrimarySpds = [];
odObsRefSpds = [];

for i = 1:15
    lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
        coneShiftData.p.Results.adjustmentLength,coneShiftData.p1,coneShiftData.p2,...
        coneShiftData.test(i),1,coneShiftData.p.Results.p2Scale,...
        coneShiftData.p.Results.testScale);
    lightData = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops',lightFile));
    
    % OD primary spds 
    sortedDiffsArr = sort(abs(lightData.p1Scales...
        -coneShiftData.primaryRatiosSim(1,i)));
    ind1 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(1,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(1,i))==sortedDiffsArr(2));
    odObsPrimarySpds = [odObsPrimarySpds,...
        mean([lightData.primarySpdsPredicted(:,ind1),lightData.primarySpdsPredicted(:,ind2)],2)];
    
    % OD ref spds 
    sortedDiffsArr = sort(abs(lightData.testScales...
        -coneShiftData.testIntensitiesSim(1,i)));
    ind1 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(1,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(1,i))==sortedDiffsArr(2));
    odObsRefSpds = [odObsRefSpds,...
        mean([lightData.testSpdsPredicted(:,ind1),lightData.testSpdsPredicted(:,ind2)],2)];
    
    % Lambda max primary spds 
    sortedDiffsArr = sort(abs(lightData.p1Scales...
        -coneShiftData.primaryRatiosSim(1,i)));
    ind1 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(2,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(2,i))==sortedDiffsArr(2));
    lambdaMaxObsPrimarySpds = [lambdaMaxObsPrimarySpds,...
        mean([lightData.primarySpdsPredicted(:,ind1),lightData.primarySpdsPredicted(:,ind2)],2)];
    
    % Lambda max ref spds 
    sortedDiffsArr = sort(abs(lightData.testScales...
        -coneShiftData.testIntensitiesSim(1,i)));
    ind1 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(2,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(2,i))==sortedDiffsArr(2));
    lambdaMaxObsRefSpds = [lambdaMaxObsRefSpds,...
        mean([lightData.testSpdsPredicted(:,ind1) lightData.testSpdsPredicted(:,ind2)],2)];
end

% Compute excitations
lambdaMaxObsRecPrimaryExcitations = lambdaMaxObsRec.T_cones*lambdaMaxObsPrimarySpds;
lambdaMaxObsRecRefExcitations = lambdaMaxObsRec.T_cones*lambdaMaxObsRefSpds;
lambdaMaxObsPrimaryExcitations = lambdaMaxObs.T_cones*lambdaMaxObsPrimarySpds;
lambdaMaxObsRefExcitations = lambdaMaxObs.T_cones*lambdaMaxObsRefSpds;
lambdaMaxObsStdPrimaryExcitations = stdObs.T_cones*lambdaMaxObsPrimarySpds;
lambdaMaxObsStdRefExcitations = stdObs.T_cones*lambdaMaxObsRefSpds;

odObsRecPrimaryExcitations = odObsRec.T_cones*odObsPrimarySpds;
odObsRecRefExcitations = odObsRec.T_cones*odObsRefSpds;
odObsPrimaryExcitations = odObs.T_cones*odObsPrimarySpds;
odObsRefExcitations = odObs.T_cones*odObsRefSpds;
odObsStdPrimaryExcitations = stdObs.T_cones*odObsPrimarySpds;
odObsStdRefExcitations = stdObs.T_cones*odObsRefSpds;
coneExcitations = {lambdaMaxObsPrimaryExcitations lambdaMaxObsRefExcitations ...
    lambdaMaxObsRecPrimaryExcitations lambdaMaxObsRecRefExcitations ...
    lambdaMaxObsStdPrimaryExcitations lambdaMaxObsStdRefExcitations ...
    odObsPrimaryExcitations odObsRefExcitations...
    odObsRecPrimaryExcitations,odObsRecRefExcitations...
    odObsStdPrimaryExcitations,odObsStdRefExcitations};

% Compute cone contrasts
LPlusM = {};
LMinusM = {};
for i = 1:12
    LPlusM{i} = coneExcitations{i}(1,:)+coneExcitations{i}(2,:);
    LMinusM{i} = coneExcitations{i}(1,:)-coneExcitations{i}(2,:);
end 

% Set up figure and subplots
coneContrastSim = figure(20);
hold on;
theAxesGrid = plotlab.axesGrid(coneContrastSim, ...
    'rowsNum', 1, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Plot results
% Figure 1: lambda max obs
set(gcf,'CurrentAxes',theAxesGrid{1,1});
axis(theAxesGrid{1,1});
xlim([0 0.6]);
ylim([0.02 0.026]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title('L Lambda Max Shift - 4 nm');
hold on;
stdPrimaryLine  = plot(LMinusM{5}./LPlusM{5},...
    LPlusM{5},'* ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
stdRefLine = plot(LMinusM{6}./LPlusM{6},...
    LPlusM{6},'o ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
simPrimaryLine  = plot(LMinusM{1}./LPlusM{1},...
    LPlusM{1},'* ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
simRefLine = plot(LMinusM{2}./LPlusM{2},...
    LPlusM{2},'o ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
recPrimaryLine  = plot(LMinusM{3}./LPlusM{3},...
    LPlusM{3},'* ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
recRefLine = plot(LMinusM{4}./LPlusM{4},...
    LPlusM{4},'o ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
for k = 1:length(LMinusM{1})
    plot([LMinusM{2}(k)./LPlusM{2}(k),LMinusM{1}(k)./LPlusM{1}(k)],...
        [LPlusM{2}(k) LPlusM{1}(k)],'k-');
     plot([LMinusM{4}(k)./LPlusM{4}(k),LMinusM{3}(k)./LPlusM{3}(k)],...
        [LPlusM{4}(k) LPlusM{3}(k)],'k-');
     plot([LMinusM{6}(k)./LPlusM{6}(k),LMinusM{5}(k)./LPlusM{5}(k)],...
        [LPlusM{6}(k) LPlusM{5}(k)],'k-');
end
legend([stdRefLine simRefLine recRefLine], 'Standard', 'Simulated', 'Recovered');

% Figure 2: OD rec
set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2});
xlim([0 0.6]);
ylim([0.02 0.026]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title('L Optical Density Shift - 18%', 'Interpreter', 'none');
hold on;
stdPrimaryLine  = plot(LMinusM{11}./LPlusM{11},...
    LPlusM{11},'* ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
stdRefLine = plot(LMinusM{12}./LPlusM{12},...
    LPlusM{12},'o ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
simPrimaryLine  = plot(LMinusM{7}./LPlusM{7},...
    LPlusM{7},'* ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
simRefLine = plot(LMinusM{8}./LPlusM{8},...
    LPlusM{8},'o ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
recPrimaryLine  = plot(LMinusM{9}./LPlusM{9},...
    LPlusM{9},'* ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
recRefLine = plot(LMinusM{10}./LPlusM{10},...
    LPlusM{10},'o ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
for k = 1:length(LMinusM{1})
    plot([LMinusM{11}(k)./LPlusM{11}(k),LMinusM{12}(k)./LPlusM{12}(k)],...
        [LPlusM{11}(k) LPlusM{12}(k)],'k-');
     plot([LMinusM{8}(k)./LPlusM{8}(k),LMinusM{7}(k)./LPlusM{7}(k)],...
        [LPlusM{8}(k) LPlusM{7}(k)],'k-');
     plot([LMinusM{10}(k)./LPlusM{10}(k),LMinusM{9}(k)./LPlusM{9}(k)],...
        [LPlusM{10}(k) LPlusM{9}(k)],'k-');
end
legend([stdRefLine simRefLine recRefLine], 'Standard', 'Simulated', 'Recovered');

sgtitle('Cone Response Difference - Noiseless Observers','FontSize',18);

%% Rayleigh Simulation Plot - Cone L and M Recovery Figure (Noisy)
% Generate observers with the desired parameters
lODParams = [0 0 2 0 0 0 0 0].*coneParamSds;
lLambdaMaxParams = [0 0 0 0 0 2 0 0].*coneParamSds;
S = [380 2 201];
opponentParams = [40.3908 205.7353 62.9590 1.0000];
odObs = genRayleighObserver('coneVec',lODParams,'S',S,'opponentParams',...
    opponentParams);
lambdaMaxObs = genRayleighObserver('coneVec',lLambdaMaxParams,'S',S,...
    'opponentParams',opponentParams);
stdObs = genRayleighObserver('coneVec',zeros(1,8),'S',S,...
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

% Find spds
lambdaMaxObsPrimarySpds = [];
lambdaMaxObsRefSpds = [];
odObsPrimarySpds = [];
odObsRefSpds = [];

for i = 1:15
    lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
        coneShiftData.p.Results.adjustmentLength,coneShiftData.p1,coneShiftData.p2,...
        coneShiftData.test(i),1,coneShiftData.p.Results.p2Scale,...
        coneShiftData.p.Results.testScale);
    lightData = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops',lightFile));
    
    % OD primary spds 
    sortedDiffsArr = sort(abs(lightData.p1Scales...
        -coneShiftData.primaryRatiosSim(1,i)));
    ind1 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(1,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(1,i))==sortedDiffsArr(2));
    odObsPrimarySpds = [odObsPrimarySpds,...
        mean([lightData.primarySpdsPredicted(:,ind1),lightData.primarySpdsPredicted(:,ind2)],2)];
    
    % OD ref spds 
    sortedDiffsArr = sort(abs(lightData.testScales...
        -coneShiftData.testIntensitiesSim(1,i)));
    ind1 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(1,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(1,i))==sortedDiffsArr(2));
    odObsRefSpds = [odObsRefSpds,...
        mean([lightData.testSpdsPredicted(:,ind1),lightData.testSpdsPredicted(:,ind2)],2)];
    
    % Lambda max primary spds 
    sortedDiffsArr = sort(abs(lightData.p1Scales...
        -coneShiftData.primaryRatiosSim(1,i)));
    ind1 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(2,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.p1Scales-coneShiftData.primaryRatiosSim(2,i))==sortedDiffsArr(2));
    lambdaMaxObsPrimarySpds = [lambdaMaxObsPrimarySpds,...
        mean([lightData.primarySpdsPredicted(:,ind1),lightData.primarySpdsPredicted(:,ind2)],2)];
    
    % Lambda max ref spds 
    sortedDiffsArr = sort(abs(lightData.testScales...
        -coneShiftData.testIntensitiesSim(1,i)));
    ind1 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(2,i))==sortedDiffsArr(1));
    ind2 = find(abs(lightData.testScales-coneShiftData.testIntensitiesSim(2,i))==sortedDiffsArr(2));
    lambdaMaxObsRefSpds = [lambdaMaxObsRefSpds,...
        mean([lightData.testSpdsPredicted(:,ind1) lightData.testSpdsPredicted(:,ind2)],2)];
end

% Compute excitations
lambdaMaxObsRecPrimaryExcitations = lambdaMaxObsRec.T_cones*lambdaMaxObsPrimarySpds;
lambdaMaxObsRecRefExcitations = lambdaMaxObsRec.T_cones*lambdaMaxObsRefSpds;
lambdaMaxObsPrimaryExcitations = lambdaMaxObs.T_cones*lambdaMaxObsPrimarySpds;
lambdaMaxObsRefExcitations = lambdaMaxObs.T_cones*lambdaMaxObsRefSpds;
lambdaMaxObsStdPrimaryExcitations = stdObs.T_cones*lambdaMaxObsPrimarySpds;
lambdaMaxObsStdRefExcitations = stdObs.T_cones*lambdaMaxObsRefSpds;

odObsRecPrimaryExcitations = odObsRec.T_cones*odObsPrimarySpds;
odObsRecRefExcitations = odObsRec.T_cones*odObsRefSpds;
odObsPrimaryExcitations = odObs.T_cones*odObsPrimarySpds;
odObsRefExcitations = odObs.T_cones*odObsRefSpds;
odObsStdPrimaryExcitations = stdObs.T_cones*odObsPrimarySpds;
odObsStdRefExcitations = stdObs.T_cones*odObsRefSpds;
coneExcitations = {lambdaMaxObsPrimaryExcitations lambdaMaxObsRefExcitations ...
    lambdaMaxObsRecPrimaryExcitations lambdaMaxObsRecRefExcitations ...
    lambdaMaxObsStdPrimaryExcitations lambdaMaxObsStdRefExcitations ...
    odObsPrimaryExcitations odObsRefExcitations...
    odObsRecPrimaryExcitations,odObsRecRefExcitations...
    odObsStdPrimaryExcitations,odObsStdRefExcitations};

% Compute cone contrasts
LPlusM = {};
LMinusM = {};
for i = 1:12
    LPlusM{i} = coneExcitations{i}(1,:)+coneExcitations{i}(2,:);
    LMinusM{i} = coneExcitations{i}(1,:)-coneExcitations{i}(2,:);
end 

% Set up figure and subplots
coneContrastSimNoisy = figure(21);
hold on;
theAxesGrid = plotlab.axesGrid(coneContrastSimNoisy, ...
    'rowsNum', 1, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Plot results
% Figure 1: lambda max obs
set(gcf,'CurrentAxes',theAxesGrid{1,1});
axis(theAxesGrid{1,1});
xlim([0 0.6]);
ylim([0.02 0.026]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title('L Lambda Max Shift - 4 nm');
hold on;
stdPrimaryLine  = plot(LMinusM{5}./LPlusM{5},...
    LPlusM{5},'* ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
stdRefLine = plot(LMinusM{6}./LPlusM{6},...
    LPlusM{6},'o ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
simPrimaryLine  = plot(LMinusM{1}./LPlusM{1},...
    LPlusM{1},'* ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
simRefLine = plot(LMinusM{2}./LPlusM{2},...
    LPlusM{2},'o ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
recPrimaryLine  = plot(LMinusM{3}./LPlusM{3},...
    LPlusM{3},'* ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
recRefLine = plot(LMinusM{4}./LPlusM{4},...
    LPlusM{4},'o ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
for k = 1:length(LMinusM{1})
    plot([LMinusM{2}(k)./LPlusM{2}(k),LMinusM{1}(k)./LPlusM{1}(k)],...
        [LPlusM{2}(k) LPlusM{1}(k)],'k-');
     plot([LMinusM{4}(k)./LPlusM{4}(k),LMinusM{3}(k)./LPlusM{3}(k)],...
        [LPlusM{4}(k) LPlusM{3}(k)],'k-');
     plot([LMinusM{6}(k)./LPlusM{6}(k),LMinusM{5}(k)./LPlusM{5}(k)],...
        [LPlusM{6}(k) LPlusM{5}(k)],'k-');
end
legend([stdRefLine simRefLine recRefLine], 'Standard', 'Simulated', 'Recovered');

% Figure 2: OD rec
set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2});
xlim([0 0.6]);
ylim([0.02 0.026]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title('L Optical Density Shift - 18%', 'Interpreter', 'none');
hold on;
stdPrimaryLine  = plot(LMinusM{11}./LPlusM{11},...
    LPlusM{11},'* ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
stdRefLine = plot(LMinusM{12}./LPlusM{12},...
    LPlusM{12},'o ','MarkerEdgeColor',...
    'Red','MarkerFaceColor','Red','MarkerSize',10);
simPrimaryLine  = plot(LMinusM{7}./LPlusM{7},...
    LPlusM{7},'* ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
simRefLine = plot(LMinusM{8}./LPlusM{8},...
    LPlusM{8},'o ','MarkerEdgeColor',...
    'Green','MarkerFaceColor','Green','MarkerSize',7);
recPrimaryLine  = plot(LMinusM{9}./LPlusM{9},...
    LPlusM{9},'* ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
recRefLine = plot(LMinusM{10}./LPlusM{10},...
    LPlusM{10},'o ','MarkerEdgeColor',...
    'Blue','MarkerFaceColor','Blue','MarkerSize',5);
for k = 1:length(LMinusM{1})
    plot([LMinusM{11}(k)./LPlusM{11}(k),LMinusM{12}(k)./LPlusM{12}(k)],...
        [LPlusM{11}(k) LPlusM{12}(k)],'k-');
     plot([LMinusM{8}(k)./LPlusM{8}(k),LMinusM{7}(k)./LPlusM{7}(k)],...
        [LPlusM{8}(k) LPlusM{7}(k)],'k-');
     plot([LMinusM{10}(k)./LPlusM{10}(k),LMinusM{9}(k)./LPlusM{9}(k)],...
        [LPlusM{10}(k) LPlusM{9}(k)],'k-');
end
legend([stdRefLine simRefLine recRefLine], 'Standard', 'Simulated', 'Recovered');

sgtitle('Cone Response Difference - Noisy Observers','FontSize',18);