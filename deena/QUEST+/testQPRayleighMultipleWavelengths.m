function testQPRayleighMultipleWavelengths(subjID,nObservers,nTrials,...
    p1Wl,p2Wl,testWls,varargin)
% Runs Quest+ Rayleigh matching experiments for test wavelengths alone or
% in combination
%
% Syntax:
%   testQPRayleighMultipleWavelengths(subjID,nObservers,nTrials,...
%   p1Wl,p2Wl,testWls)
%
% Description:
%    Tests whether an enhanced QUEST+ Rayleigh matching routine with 
%    multiple test wavelengths can recover observer parameters better than 
%    Rayleigh matches with just one test wavelength. Given a series of test
%    wavelengths, runs a QUEST+ Rayleigh matching simulation with all of
%    the test wavelengths available, as well as with each test wavelength
%    individually. Produces a plot to compare parameter recovery for the
%    different simulation approaches. 
%
% Inputs:
%    subjID            -Character vector of subject ID
%    nObservers        -Number of simulated observers to test.
%    nTrials           -Number of QUEST+ trials to run in each session.
%    p1Wl              -Integer value of desired wavelength for the first
%                       primary.
%    p2Wl              -Integer value of desired wavelength for the second
%                       primary.
%    testWls           -Integer or numeric vector of desired wavelengths
%                       for the test light.
% Outputs: 
%    None (produces a figure)
%
% Optional key-value pairs:
%    'fieldSize'         -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'opponentParams'    -4-element vector with opponent contrast
%                         parameters. (1) is the luminance weight, (2) is
%                         the RG weight, (3) is the BY weight, and (4) is
%                         the baseline noise standard deviation. Default is
%                         [40.3908  205.7353   62.9590    1.0000].
%    'baseConeParams'    -Eight-element numeric vector of individual
%                         difference parameters: lens pigment density,
%                         macular pigment density, L/M/S photopigment
%                         densities, and L/M/S peak spectral sensitivities
%                         (lambda max). The values are used as means for
%                         observer sampling. Default is zeros(1,8).
%    'coneParamsToVary'  -Eight-element numeric vector of ones and zeros
%                         indicating which individual difference parameters
%                         should be varied (the noise parameter is excluded).
%                         Parameters set to 1 will be sampled around their
%                         means, while parameters set to 0 will stay at the
%                         values specified in baseParams. Default is 
%                         [0 0 1 1 0 1 1 0].
%    'noiseScaleFactor'  -Number >=0 which determines which scalar multiple
%                         of the opponent noise SD should be used as the
%                         observer noise SD. Default is 0.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'nStimValues'       -Number of possible stimulus values for lambda and
%                         test intensity, spaced evenly between 0 and 1.
%                         Default is 301.
%    'nPsiValues'        -Number of possible levels for the recovered cone
%                         parameters, spaced evenly between -2 and 2
%                         standard deviations. Default is 8.
%    'lambdaRef'         -Numerical scale factor between 0 and 1 for the
%                         value of lambda used for the reference light,
%                         which is used as a baseline for computing
%                         opponent contrasts. If set to [], no reference 
%                         light is used, and the opponent contrast of the
%                         test is calculated relative to the primary 
%                         mixture. Default is [].
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on
%                         stimulus parameters. Each row represents a given
%                         test wavelength and the limits which are associated
%                         with it. The columns are arranged as follows:
%                         [test wl, min lambda, max lambda, min test
%                         intensity, max test intensity]. Default is 
%                         [570.0000   0    0.1264    0.0399    0.0459;...
%                          590.0000    0.0456    0.4746    0.0462    0.0716;...
%                          610.0000    0.2617    0.8120    0.0695    0.1325;...
%                          630.0000    0.6046    0.9604    0.1619    0.2685;...
%                          650.0000    0.8688    0.9938    0.5109    0.6458]

% History
%    12/21/20   dce   - Wrote it
%    12/22/20   dce   - Added file saving
%    12/28/20   dce   - Separated L and M cone variation
%    12/30/20   dce   - Fixed plotting
%
% Example:
%    testQPRayleighMultipleWavelengths('test1vsMany',1,170,670,560,[570 590
%    610 630 650])
%

% Parse input
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('baseConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('coneParamsToVary',[0 0 1 1 0 1 1 0],@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('noiseScaleFactor',0,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('nStimValues',301,@(x)(isnumeric(x)));
p.addParameter('nPsiValues',9,@(x)(isnumeric(x)));
p.addParameter('lambdaRef',[],@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('stimLimits',[570.0000   0    0.1264    0.0399    0.0459;...
    590.0000    0.0456    0.4746    0.0462    0.0716;...
    610.0000    0.2617    0.8120    0.0695    0.1325;...
    630.0000    0.6046    0.9604    0.1619    0.2685;...
    650.0000    0.8688    0.9938    0.5109    0.6458],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Create output directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'Quest',subjID);
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end 

% Sample observers 
sampledConeParamsL = sampleRayleighObservers(nObservers,p.Results.baseConeParams,...
    [0 0 1 0 0 1 0 0]);
sampledConeParamsM = sampleRayleighObservers(nObservers,p.Results.baseConeParams,...
    [0 0 0 1 0 0 1 0]);

% Test QUEST+ using all test wavelengths
[~,~,recoveredParamsCombL] = qpRayleighSim([subjID '_L'],1,nObservers,nTrials,...
    p.Results.baseConeParams,[0 0 1 0 0 1 0 0],...
    p.Results.noiseScaleFactor,p1Wl,p2Wl,testWls,'precomputeQuest',false,...
    'age',p.Results.age,'fieldSize',p.Results.fieldSize,'stimLimits',...
    p.Results.stimLimits,'sampledObservers',sampledConeParamsL,'S',...
    p.Results.S,'lambdaRef',p.Results.lambdaRef,'nStimValues',...
    p.Results.nStimValues,'nPsiValues',p.Results.nPsiValues,'opponentParams',...
    p.Results.opponentParams,'p1Scale',p.Results.p1Scale,'p2Scale',...
    p.Results.p2Scale,'testScale',p.Results.testScale,'plotAll',false,...
    'plotLast',false);
[~,~,recoveredParamsCombM] = qpRayleighSim([subjID '_M'],1,nObservers,nTrials,...
    p.Results.baseConeParams,[0 0 0 1 0 0 1 0],...
    p.Results.noiseScaleFactor,p1Wl,p2Wl,testWls,'precomputeQuest',false,...
    'age',p.Results.age,'fieldSize',p.Results.fieldSize,'stimLimits',...
    p.Results.stimLimits,'sampledObservers',sampledConeParamsM,'S',...
    p.Results.S,'lambdaRef',p.Results.lambdaRef,'nStimValues',...
    p.Results.nStimValues,'nPsiValues',p.Results.nPsiValues,'opponentParams',...
    p.Results.opponentParams,'p1Scale',p.Results.p1Scale,'p2Scale',...
    p.Results.p2Scale,'testScale',p.Results.testScale,'plotAll',false,...
    'plotLast',false);

% Test QUEST+ using individual test wavelengths
recoveredParamsIndL = cell(1,length(testWls));
recoveredParamsIndM = cell(1,length(testWls));
for i = 1:length(testWls)
    limMatrix = p.Results.stimLimits(p.Results.stimLimits(:,1)==testWls(i),:);
    [~,~,recoveredParamsIndL{i}] = qpRayleighSim([subjID '_L'],i+1,...
        nObservers,nTrials,p.Results.baseConeParams,[0 0 1 0 0 1 0 0],...
        p.Results.noiseScaleFactor,p1Wl,p2Wl,testWls(i),'precomputeQuest',false,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,'stimLimits',...
        limMatrix,'sampledObservers',sampledConeParamsL,'S',...
        p.Results.S,'lambdaRef',p.Results.lambdaRef,'nStimValues',...
        p.Results.nStimValues,'nPsiValues',p.Results.nPsiValues,'opponentParams',...
        p.Results.opponentParams,'p1Scale',p.Results.p1Scale,'p2Scale',...
        p.Results.p2Scale,'testScale',p.Results.testScale,'plotAll',false,...
        'plotLast',false);
    [~,~,recoveredParamsIndM{i}] = qpRayleighSim([subjID '_M'],i+1,...
        nObservers,nTrials,p.Results.baseConeParams,[0 0 0 1 0 0 1 0],...
        p.Results.noiseScaleFactor,p1Wl,p2Wl,testWls(i),'precomputeQuest',false,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,'stimLimits',...
        limMatrix,'sampledObservers',sampledConeParamsM,'S',...
        p.Results.S,'lambdaRef',p.Results.lambdaRef,'nStimValues',...
        p.Results.nStimValues,'nPsiValues',p.Results.nPsiValues,'opponentParams',...
        p.Results.opponentParams,'p1Scale',p.Results.p1Scale,'p2Scale',...
        p.Results.p2Scale,'testScale',p.Results.testScale,'plotAll',false,...
        'plotLast',false);
end

%% Parameter recovery across conditions
combParamErrL = zeros(1,nObservers);
combParamErrM = zeros(1,nObservers);
stdParamErrL = zeros(1,nObservers);
stdParamErrM = zeros(1,nObservers);
indParamErrL = zeros(length(testWls),nObservers);
indParamErrM = zeros(length(testWls),nObservers);
for i = 1:nObservers
    combParamErrL(i) = findParamRecoveryError(sampledConeParamsL(i,:),...
        recoveredParamsCombL(i,:));
    combParamErrM(i) = findParamRecoveryError(sampledConeParamsM(i,:),...
        recoveredParamsCombM(i,:));
    stdParamErrL(i) = findParamRecoveryError(sampledConeParamsL(i,:),...
        zeros(1,8));
    stdParamErrM(i) = findParamRecoveryError(sampledConeParamsM(i,:),...
        zeros(1,8));
    % Loop through test wavelengths 
    for j = 1:length(testWls)
        indParamErrL(j,i) = findParamRecoveryError(sampledConeParamsL(i,:),...
            recoveredParamsIndL{j}(i,:));
        indParamErrM(j,i) = findParamRecoveryError(sampledConeParamsM(i,:),...
            recoveredParamsIndM{j}(i,:));
    end         
end
avgParamCombErrL = mean(combParamErrL);
avgParamCombErrM = mean(combParamErrM);
avgParamStdErrL = mean(stdParamErrL);
avgParamStdErrM = mean(stdParamErrM);
avgParamIndErrL = mean(indParamErrL,2);
avgParamIndErrM = mean(indParamErrM,2);
[minParamIndErrL,minParamIndErrLInd] = min(avgParamIndErrL);
[minParamIndErrM,minParamIndErrMInd] = min(avgParamIndErrM);
bestParamWlL = testWls(minParamIndErrLInd);
bestParamWlM = testWls(minParamIndErrMInd);


%% Cone recovery across conditions
combConeErrL = zeros(1,nObservers);
combConeErrM = zeros(1,nObservers);
stdConeErrL = zeros(1,nObservers);
stdConeErrM = zeros(1,nObservers);
indConeErrL = zeros(length(testWls),nObservers);
indConeErrM = zeros(length(testWls),nObservers);
for i = 1:nObservers
    [~,combConeErrL(i),~,~] = ...
        findConeSensitivityError(sampledConeParamsL(i,:),recoveredParamsCombL(i,:),...
        'S',p.Results.S,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
        'opponentParams',p.Results.opponentParams);
    [~,~,combConeErrM(i),~] = ...
        findConeSensitivityError(sampledConeParamsM(i,:),recoveredParamsCombM(i,:),...
        'S',p.Results.S,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
        'opponentParams',p.Results.opponentParams);
    
    [~,stdConeErrL(i),~,~] = ...
        findConeSensitivityError(sampledConeParamsL(i,:),zeros(1,8),...
        'S',p.Results.S,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
        'opponentParams',p.Results.opponentParams);
    [~,~,stdConeErrM(i),~] = ...
        findConeSensitivityError(sampledConeParamsM(i,:),zeros(1,8),...
        'S',p.Results.S,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
        'opponentParams',p.Results.opponentParams);
    
    % Loop through test wavelengths 
    for j = 1:length(testWls)
        [~,indConeErrL(j,i),~,~] = ...
            findConeSensitivityError(sampledConeParamsL(i,:),recoveredParamsIndL{j}(i,:),...
            'S',p.Results.S,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
            'opponentParams',p.Results.opponentParams);
        [~,~,indConeErrM(j,i),~] = ...
            findConeSensitivityError(sampledConeParamsM(i,:),recoveredParamsIndM{j}(i,:),...
            'S',p.Results.S,'fieldSize',p.Results.fieldSize,'age',p.Results.age,...
            'opponentParams',p.Results.opponentParams);
    end 
end
avgConeCombErrL = mean(combConeErrL);
avgConeCombErrM = mean(combConeErrM);
avgConeStdErrL = mean(stdConeErrL);
avgConeStdErrM = mean(stdConeErrM);
avgConeIndErrL = mean(indConeErrL,2);
avgConeIndErrM = mean(indConeErrM,2);
[minConeIndErrL,minConeIndErrLInd] = min(avgConeIndErrL);
[minConeIndErrM,minConeIndErrMInd] = min(avgConeIndErrM);
bestConeWlL = testWls(minConeIndErrLInd);
bestConeWlM = testWls(minConeIndErrMInd);

save(fullfile(outputDir,'test1vsMany.mat'),'subjID','nObservers','nTrials',...
    'p1Wl','p2Wl','testWls','p','sampledConeParamsL','sampledConeParamsM',...
    'recoveredParamsCombL','recoveredParamsCombM','recoveredParamsIndL',...
    'recoveredParamsIndM','combConeErrL','combConeErrM','stdConeErrL',...
    'stdConeErrM','indConeErrL','indConeErrM','avgConeCombErrL',...
    'avgConeCombErrM','avgConeStdErrL','avgConeStdErrM','avgConeIndErrL',...
    'avgConeIndErrM','minConeIndErrL','minConeIndErrM','bestConeWlL',...
    'bestConeWlM','combParamErrL','combParamErrM','stdParamErrL',...
    'stdParamErrM','indParamErrL','indParamErrM','avgParamCombErrL',...
    'avgParamCombErrM','avgParamStdErrL','avgParamStdErrM','avgParamIndErrL',...
    'avgParamIndErrM','minParamIndErrL','minParamIndErrM','bestParamWlL',...
    'bestParamWlM');
%% Plots
%
% 1 - parameter recovery (all wavelengths)
% Initial setup 
nCols = 2;
nRows = 2;
coneParamNames = {'L Photopigment Density','M Photopigment Density',...
    'L Lambda Max','M Lambda Max'};
paramIndsToPlot = [3 4 6 7];
colorString = 'bgrymcbgrymc';
legendCell = cellstr(num2str(testWls'))';
legendCell(length(testWls)+1) = {'Combined'}; 

% Set up figure
paramsPlot = figure(); clf;
set(paramsPlot,'Color',[1 1 1],'Position',[10 10 800 800]);
hold on;
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', nRows, ...
    'colsNum', nCols, ...
    'heightMargin',  0.1, ...
    'widthMargin',    0.07, ...
    'leftMargin',     0.04, ...
    'rightMargin',    0.04, ...
    'bottomMargin',   0.07, ...
    'topMargin',      0.1);

% Plot results
for i = 1:nRows*nCols
    % Make a subplot in the correct position
    row = ceil(i/nCols);
    col = mod(i,nCols);
    if col == 0
        col = nCols;
    end
    subplot('Position', subplotPosVectors(row,col).v);
    hold on;
    
    % Define axis limits and which arrays we are searching in
    if (i==3) || (i==4)  % Lambda max shifts, in nm
        limits = [-5 5];
    else                   % Density shifts, in percent
        limits = [-40 40];
    end
    if (i==1) || (i==3)    % L cone
        sampledConeParams = sampledConeParamsL;
        recoveredParamsComb = recoveredParamsCombL;
        recoveredParamsInd = recoveredParamsIndL;
    else                    % M cone
        sampledConeParams = sampledConeParamsM;
        recoveredParamsComb = recoveredParamsCombM;
        recoveredParamsInd = recoveredParamsIndM;
    end 
    xlim(limits);
    ylim(limits);
    axis('square');
    
    % Plot data
    xVals = sampledConeParams(:,paramIndsToPlot(i));  % Predicted parameters
    for j = 1:length(testWls)
        yVals = recoveredParamsInd{j}(:,paramIndsToPlot(i));
        plot(xVals,yVals,'o ','MarkerSize',5,'MarkerFaceColor',...
            colorString(j),'MarkerEdgeColor',colorString(j));
    end 
    combVals = recoveredParamsComb(:,paramIndsToPlot(i)); % Recovered params
    plot(xVals,combVals,'k* ','MarkerSize',7);
    refline(1,0); 
    
    % Titles and labels
    theTitle = sprintf('%s Recovered vs Simulated',cell2mat(coneParamNames(i)));
    title(theTitle);
    xlabel('Simulated Parameters');
    ylabel('Recovered Parameters');
    legend(legendCell);
end
sgtitle('Simulated vs Recovered Parameters');
NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_'....
     '_paramRecoveryComparison.pdf']),paramsPlot,300);
 
% 2 - parameter recovery (best wavelengths) 
% Initial setup 
nCols = 2;
nRows = 2;
coneParamNames = {'L Photopigment Density','M Photopigment Density',...
    'L Lambda Max','M Lambda Max'};
paramIndsToPlot = [3 4 6 7];

% Set up figure
paramsPlotBest = figure(); clf;
set(paramsPlotBest,'Color',[1 1 1],'Position',[10 10 800 800]);
hold on;
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', nRows, ...
    'colsNum', nCols, ...
    'heightMargin',  0.1, ...
    'widthMargin',    0.07, ...
    'leftMargin',     0.04, ...
    'rightMargin',    0.04, ...
    'bottomMargin',   0.07, ...
    'topMargin',      0.1);

% Plot results
for i = 1:nRows*nCols
    % Make a subplot in the correct position
    row = ceil(i/nCols);
    col = mod(i,nCols);
    if col == 0
        col = nCols;
    end
    subplot('Position', subplotPosVectors(row,col).v);
    hold on;
    
    % Define axis limits and which arrays we are searching in
    if (i==3) || (i==4)  % Lambda max shifts, in nm
        limits = [-5 5];
    else                   % Density shifts, in percent
        limits = [-40 40];
    end
    if (i==1) || (i==3)    % L cone
        sampledConeParams = sampledConeParamsL;
        recoveredParamsComb = recoveredParamsCombL;
        recoveredParamsInd = recoveredParamsIndL;
        bestWl = bestConeWlL;
    else                    % M cone
        sampledConeParams = sampledConeParamsM;
        recoveredParamsComb = recoveredParamsCombM;
        recoveredParamsInd = recoveredParamsIndM;
        bestWl = bestConeWlM; 
    end 
    xlim(limits);
    ylim(limits);
    axis('square');
    
    % Plot data
    xVals = sampledConeParams(:,paramIndsToPlot(i));  % Predicted parameters
    for j = 1:length(testWls)
        if testWls(j)==bestWl
            yVals = recoveredParamsInd{j}(:,paramIndsToPlot(i));
            plot(xVals,yVals,'o ','MarkerSize',5,'MarkerFaceColor',...
                'Blue','MarkerEdgeColor','Blue');
        end
    end
    combVals = recoveredParamsComb(:,paramIndsToPlot(i)); % Recovered params
    plot(xVals,combVals,'k* ','MarkerSize',7);
    refline(1,0); 
    
    % Titles and labels
    theTitle = sprintf('%s Recovered vs Simulated',cell2mat(coneParamNames(i)));
    title(theTitle);
    xlabel('Simulated Parameters');
    ylabel('Recovered Parameters');
    legend(num2str(bestWl),'Combined');
end
sgtitle('Simulated vs Recovered Parameters - Best Individual Test Wavelengths');
NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_'....
     '_paramRecoveryBest.pdf']),paramsPlotBest,300); 

 
%% 3 - error plot
coneErrPlot = figure(); clf;
subplot(2,1,1);
dataArr = [avgConeCombErrL minConeIndErrL avgConeStdErrL];
bar(dataArr);
name = {'Combined','Best Individual','Standard'};
set(gca,'xticklabel',name)
title('Average L Cone Recovery Error');
ylabel('RMS Error');
ylim([0 0.02]);
for i=1:length(dataArr)
    text(i,dataArr(i),num2str(dataArr(i),'%0.4f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
end 

subplot(2,1,2);
dataArr = [avgConeCombErrM minConeIndErrM avgConeStdErrM];
bar(dataArr);
name = {'Combined','Best Individual','Standard'};
set(gca,'xticklabel',name)
title('Average M Cone Recovery Error');
ylabel('RMS Error');
ylim([0 0.02]);
for i=1:length(dataArr)
    text(i,dataArr(i),num2str(dataArr(i),'%0.4f'),...
        'HorizontalAlignment','center',...
        'VerticalAlignment','bottom');
end 

NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_'....
     '_coneErrPlot.pdf']),coneErrPlot,300); 
end