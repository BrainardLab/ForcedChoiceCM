% Generates figures for Rayleigh matching manuscript

% History 
%     11/5/20  dce   -Wrote it 
%     11/8/20  dce   -Added opponent sphere

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
    'figureHeightInches', 7,...
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
coneShiftPlot = figure(1);
hold on;
theAxesGrid = plotlab.axesGrid(coneShiftPlot, ...
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

%% Opponent contrast sphere figure figure
% May want to fix random number generator at some point
noiseSD = 1;
[~,~,opponentSphere] = opponentAxesToLab(noiseSD);
figure(3);
plot3(opponentSphere(1,:),opponentSphere(2,:),opponentSphere(3,:),'o-');
xlabel('x');
ylabel('y');
zlabel('z');
title('Optimized Opponent Sphere, Radius = 1');

%% Rayleigh simulation plots