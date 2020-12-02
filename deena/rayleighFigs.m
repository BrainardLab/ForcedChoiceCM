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
% Produces three figures. The first is a two panel figure where the first 
% panel show how match patterns change with lambda max, and the second 
% shows how this is obscured by optical density shifts. The second figure 
% shows how the patches can be separated by using different test wavelengths

% Initialize figure
% Set up figure and subplots
pittPlot = figure(2);
hold on;
theAxesGrid = plotlab.axesGrid(pittPlot, ...
    'rowsNum', 1, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

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

set(gcf,'CurrentAxes',theAxesGrid{1,2});
theTitle = 'Vary L $\lambda_{max}$ and Optical Density'; 
plotRayleighMatchesObserver(lambdaMaxODObs1,p1Spd,p2Spd,testSpds(:,testWls==testWlInitial),...
    noiseScalar,colors{1},theTitle,'figHandle',pittPlot)
plotRayleighMatchesObserver(lambdaMaxODObs2,p1Spd,p2Spd,testSpds(:,testWls==testWlInitial),...
    noiseScalar,colors{2},theTitle,'figHandle',pittPlot)
legend('$\lambda_{max}$ +2 Sd, OD -2 Sd','OD +2 Sd','interpreter','latex')

gPittPlot = figure(3);
theTitle = 'Generalized Pitt Diagram: Vary L $\lambda_{max}$ and Optical Density'; 
for i = 1:length(testWls)
    plotRayleighMatchesObserver(lambdaMaxODObs1,p1Spd,p2Spd,testSpds(:,i),...
        noiseScalar,colors{1},theTitle,'figHandle',gPittPlot)
    plotRayleighMatchesObserver(lambdaMaxODObs2,p1Spd,p2Spd,testSpds(:,i),...
        noiseScalar,colors{2},theTitle,'figHandle',gPittPlot)
end 
legend('$\lambda_{max}$ +2 Sd, OD -2 Sd','OD +2 Sd','interpreter','latex')

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
legend('Test','Location', 'NorthWest');
title('Test Spd (600 nm)'); 
xlabel('Wavelength (nm)');
ylabel('Radiance (W/[sr-m^2-nm])');

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
[~,~,~,~,~,~,recoveredParams]...
    = sampleRayleighMatch('coneShift',2,zeros(1,8),[0 0 1 1 0 1 1 0],...
    opponentParams,670,540,570:5:640,'forcedChoice','noiseScaleFactor',...
    0,'sampledObservers',[lODParams;lLambdaMaxParams],'makeNoObserverPlots',...
    true);
odObsRec = genRayleighObserver('coneVec',recoveredParams(1,:),'S',S,...
    'opponentParams',opponentParams);
lambdaMaxObsRec = genRayleighObserver('coneVec',recoveredParams(2,:),'S',S,...
    'opponentParams',opponentParams);

% Set up figure and subplots
coneRecoveryPlotNoiseless = figure(6);
hold on;
theAxesGrid = plotlab.axesGrid(coneRecoveryPlotNoiseless, ...
    'rowsNum', 3, 'colsNum', 2, ...
    'spacing', 'normal', ...
    'padding', 'normal', ...
    'method', 'tile');

% Plot results
set(gcf,'CurrentAxes',theAxesGrid{1,1});
axis(theAxesGrid{1,1}, 'square');
plot(wls,standardObs.T_cones(1,:),'g',wls,lambdaMaxObs.T_cones(1,:),'b',...
    wls,lambdaMaxObsRec.T_cones(1,:));
theTitle = sprintf('L cone Lambda Max Shift: %g nm',coneParamSds(6)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated', 'Recovered', 'FontSize',12);
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.015, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{1,2});
axis(theAxesGrid{1,2}, 'square');
plot(wls,standardObs.T_cones(1,:),'g',wls,odObs.T_cones(1,:),'b',...
    wls,odObsRec.T_cones(1,:),'r');
theTitle = sprintf('L cone Optical Density Shift: %g%%',coneParamSds(3)*2);
title(theTitle);
lgd = legend('Standard', 'Simulated', 'Recovered', 'FontSize',12);
plotlab.repositionLegend(lgd, [lgd.Position(1)+0.02, lgd.Position(2)])
xlabel('Wavelength (nm)');
ylabel('Sensitivity');

set(gcf,'CurrentAxes',theAxesGrid{2,1});
axis(theAxesGrid{2,1}, 'square');
plot(wls,lambdaMaxObs.T_cones(1,:)-standardObs.T_cones(1,:),'g',...
    wls,lambdaMaxObs.T_cones(1,:)-lambdaMaxObsRec.T_cones(1,:),'r');
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Cone Sensitivity Difference');
ylim([-0.06 0.06]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{2,2}, 'square');
plot(wls,odObs.T_cones(1,:)-standardObs.T_cones(1,:),'g',...
    wls,odObs.T_cones(1,:)-odObsRec.T_cones(1,:),'r');
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Cone Sensitivity Difference');

set(gcf,'CurrentAxes',theAxesGrid{3,1});
axis(theAxesGrid{2,1}, 'square');
plot(wls,lambdaMaxObs.T_cones(2,:)-standardObs.T_cones(2,:),'g',...
    wls,lambdaMaxObs.T_cones(2,:)-lambdaMaxObsRec.T_cones(2,:),'r');
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('M Cone Sensitivity Difference');
ylim([-0.06 0.06]);

set(gcf,'CurrentAxes',theAxesGrid{2,2});
axis(theAxesGrid{3,2}, 'square');
plot(wls,odObs.T_cones(2,:)-standardObs.T_cones(2,:),'g',...
    wls,odObs.T_cones(2,:)-odObsRec.T_cones(2,:),'r');
xlabel('Wavelength (nm)');
ylabel('Sensitivity Difference');
title('L Cone Sensitivity Difference');
ylim([-0.06 0.06]);
