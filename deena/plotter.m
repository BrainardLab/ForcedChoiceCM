% Variable params plots 
paramsNoiseFig = figure();
set(paramsNoiseFig,'Color',[1 1 1],'Position',[10 10 800 800]);
hold on;
nCols2 = 3;
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
plotNames = {'Noise SD = 0','Noise SD = 0.02',...
    'Noise SD = 0.04','Noise SD = 0.06','Noise SD = 0.08'};
file0 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNoiseForPaper2_0_FC\varyNoiseForPaper2_0_FC_paramsSearchData.mat');
file1 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNoiseForPaper2_1_FC\varyNoiseForPaper2_1_FC_paramsSearchData.mat');
file2 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNoiseForPaper2_2_FC\varyNoiseForPaper2_2_FC_paramsSearchData.mat');
file3 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNoiseForPaper2_3_FC\varyNoiseForPaper2_3_FC_paramsSearchData.mat');
file4 = load('\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\paramsSearch\varyNoiseForPaper2_4_FC\varyNoiseForPaper2_4_FC_paramsSearchData.mat');
files = {file0 file1 file2 file3 file4}; 

for i = 1:5
    % Make a subplot in the correct position
    row = ceil(i/nCols2);
    col = mod(i,nCols2);
    if col == 0
        col = nCols2;
    end
    subplot('Position', subplotPosVectors(row,col).v);
    hold on;
    xlim([-5 5]);
    ylim([-5 5]);
    axis('square');
    theFile = cell2mat(files(i));
    xVals = theFile.sampledParams(:,7);   % Predicted parameters
    yVals = theFile.recoveredParams(:,7);         % Recovered params
    l1 = plot(xVals,yVals,'b* ','MarkerSize',7,'LineWidth',1);
    l2 = refline(1,0);
    theTitle = sprintf(cell2mat(plotNames(i)));
    title(theTitle);
    xlabel('Simulated Parameter');
    ylabel('Recovered Parameter');
    lgd = legend([l1 l2],'Parameter','y=x');
    lgd.Location = 'northwest';
end 
sgtitle('M Cone Lambda Max Recovery - Vary Noise')

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