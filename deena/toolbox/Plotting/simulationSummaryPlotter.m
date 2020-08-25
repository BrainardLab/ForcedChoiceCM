% Variable params plots 
paramsNoiseFig = figure();
set(paramsNoiseFig,'Color',[1 1 1],'Position',[10 10 800 800]);
hold on;
nCols2 = 2;
nRows2 = 2;
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', nRows2, ...
    'colsNum', nCols2, ...
    'heightMargin',  0.1, ...
    'widthMargin',    0.07, ...
    'leftMargin',     0.04, ...
    'rightMargin',    0.04, ...
    'bottomMargin',   0.1, ...
    'topMargin',      0.1);
plotNames = {'1 Match','2 Matches',...
    '3 Matches','4 Matches'};
file1 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNMatchesForPaper7_1_FC\varyNMatchesForPaper7_1_FC_paramsSearchData.mat');
file2 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNMatchesForPaper7_2_FC\varyNMatchesForPaper7_2_FC_paramsSearchData.mat');
file3 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNMatchesForPaper7_3_FC\varyNMatchesForPaper7_3_FC_paramsSearchData.mat');
file4 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNMatchesForPaper7_4_FC\varyNMatchesForPaper7_4_FC_paramsSearchData.mat');
files = {file1 file2 file3 file4}; 

for i = 1:4
    % Make a subplot in the correct position
    row = ceil(i/nCols2);
    col = mod(i,nCols2);
    if col == 0
        col = nCols2;
    end
    subplot('Position', subplotPosVectors(row,col).v);
    hold on;
    xlim([-40 40]);
    ylim([-40 40]);
    axis('square');
    theFile = cell2mat(files(i));
    xVals = theFile.sampledParams(:,4);   % Predicted parameters
    yVals = theFile.recoveredParams(:,4);         % Recovered params
    l1 = plot(xVals,yVals,'b* ','MarkerSize',7,'LineWidth',1);
    l2 = refline(1,0);
    theTitle = sprintf(cell2mat(plotNames(i)));
    title(theTitle);
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
end 
sgtitle('M Cone Optical Density Recovery - Effect of Number of Matches')

%% Error chart
    plottingInds = [bestObs worstObs];
    nRows = 2;
    nCols = 2;

% Set up the base figure
errFig = figure();
set(errFig,'Color',[1 1 1],'Position',[10 10 1400 800]);
hold on;
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', nRows, ...
    'colsNum', nCols, ...
    'heightMargin',  0.1, ...
    'widthMargin',    0.1, ...
    'leftMargin',     0.07, ...
    'rightMargin',    0.04, ...
    'bottomMargin',   0.1, ...
    'topMargin',      0.1);

% Subplot 1 - cone spectral sensitivity error
subplot('Position', subplotPosVectors(1,1).v);
hold on;
bar([coneErr coneStandardErr]);
title('Cone Spectral Sensitivity Error');
legend('Simulated vs Recovered Params', 'Simulated vs Standard Params');
ylabel('Error');
xlabel('Observer');
ylim([0 0.02]);

% Subplot 2 - match error
subplot('Position', subplotPosVectors(2,1).v);
hold on;
bar([matchSampledErr matchErr matchStandardErr]);
title('Match Spectral Sensitivity Error');
legend('Simulated Parameters','Recovered Parameters','Standard Parameters');
ylabel('Error');
xlabel('Observer');
ylim([0 0.2]);

% Subplots 3-4: best and worst cones 
S = [380 2 201];
for k = 1:2
    subplot('Position', subplotPosVectors(k,2).v);
    hold on;
    % Find and plot observer cone fundamentals
    sampledObserver = genRayleighObserver('age',p.Results.age,...
        'fieldSize',p.Results.fieldSize,'coneVec',...
        sampledParams(plottingInds(k),:),'opponentParams',opponentParams);
    recoveredObserver = genRayleighObserver('age',p.Results.age,...
        'fieldSize',p.Results.fieldSize,'coneVec',...
        recoveredParams(plottingInds(k),:),'opponentParams',opponentParams);
    l1 = plot(SToWls(S),sampledObserver.T_cones(1:2,:),'b-',...
        'LineWidth',2.5);
    l2 = plot(SToWls(S),recoveredObserver.T_cones(1:2,:),'r-',...
        'LineWidth',1.25);
    % Clean up plot
    legend([l1(1) l2(1)],'Simulated Observer','Recovered Observer');
    xlabel('Wavelength (nm)');
    ylabel('Power');
end
% Edit titles
    subplot('Position', subplotPosVectors(1,2).v);
    title('L and M Cone Spectral Sensitivities, Best Observer');
    
    subplot('Position', subplotPosVectors(2,2).v);
    title('L and M Cone Spectral Sensitivities, Worst Observer');