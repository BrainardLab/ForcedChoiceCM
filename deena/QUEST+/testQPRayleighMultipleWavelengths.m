function testQPRayleighMultipleWavelengths(subjID,nObservers,nTrials,...
    p1Wl,p2Wl,testWls,varargin)
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
%                         observer sampling.
%    'coneParamsToVary'  -Eight-element numeric vector of ones and zeros
%                         indicating which individual difference parameters
%                         should be varied (the noise parameter is excluded).
%                         Parameters set to 1 will be sampled around their
%                         means, while parameters set to 0 will stay at the
%                         values specified in baseParams.
%    'noiseScaleFactor'  -Number >=0 which determines which scalar multiple
%                         of the opponent noise SD should be used as the
%                         observer noise SD.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'nStimValues'       -Number of possible stimulus values for lambda and
%                         test intensity, spaced evenly between 0 and 1.
%                         Default is 21.
%    'nPsiValues'        -Number of possible levels for the recovered cone
%                         parameters, spaced evenly between -2 and 2
%                         standard deviations. Default is 8.
%    'lambdaRef'         -Numerical scale factor between 0 and 1 for the
%                         value of lambda used for the reference light,
%                         which is used as a baseline for computing
%                         opponent contrasts. Default is 0.8. If set to [],
%                         no reference light is used, and the opponent
%                         contrast of the test is calculated relative to
%                         the primary mixture.
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on
%                         stimulus parameters. Each row represents a given
%                         test wavelength and the limits which are associated
%                         with it. The columns are arranged as follows:
%                         [test wl, min lambda, max lambda, min test
%                         intensity, max test intensity]. Default is [].
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
p.addParameter('nPsiValues',8,@(x)(isnumeric(x)));
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
sampledConeParams = sampleRayleighObservers(nObservers,p.Results.baseConeParams,...
    p.Results.coneParamsToVary);

% Test QUEST+ using all test wavelengths
[~,~,recoveredParamsComb] = qpRayleighSim(subjID,1,nObservers,nTrials,...
    p.Results.baseConeParams,p.Results.coneParamsToVary,...
    p.Results.noiseScaleFactor,p1Wl,p2Wl,testWls,'precomputeQuest',true,...
    'age',p.Results.age,'fieldSize',p.Results.fieldSize,'stimLimits',...
    p.Results.stimLimits,'sampledObservers',sampledConeParams,'S',...
    p.Results.S,'lambdaRef',p.Results.lambdaRef,'nStimValues',...
    p.Results.nStimValues,'nPsiValues',p.Results.nPsiValues,'opponentParams',...
    p.Results.opponentParams,'p1Scale',p.Results.p1Scale,'p2Scale',...
    p.Results.p2Scale,'testScale',p.Results.testScale,'plotAll',false,...
    'plotLast',false);

% Test QUEST+ using individual test wavelengths
recoveredParamsInd = cell(1,length(testWls));
for i = 1:length(testWls)
    [~,recoveredParamsInd{i}] = qpRayleighSim(subjID,i+1,nObservers,nTrials,...
        p.Results.baseConeParams,p.Results.coneParamsToVary,...
        p.Results.noiseScaleFactor,p1Wl,p2Wl,testWls(i),'precomputeQuest',false,...
        'age',p.Results.age,'fieldSize',p.Results.fieldSize,'stimLimits',...
        p.Results.stimLimits,'sampledObservers',sampledConeParams,'S',...
        p.Results.S,'lambdaRef',p.Results.lambdaRef,'nStimValues',...
        p.Results.nStimValues,'nPsiValues',p.Results.nPsiValues,'opponentParams',...
        p.Results.opponentParams,'p1Scale',p.Results.p1Scale,'p2Scale',...
        p.Results.p2Scale,'testScale',p.Results.testScale,'plotAll',false,...
        'plotLast',false);
end

%% Plot parameter recovery across conditions
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
set(paramsPlot,'Color',[1 1 1],'Position',[10 10 1700 800]);
hold on;
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
    'rowsNum', nRows, ...
    'colsNum', nCols, ...
    'heightMargin',  0.07, ...
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
    
    % Define axis limits
    if (i==3) || (i==4)  % Lambda max shifts, in nm
        limits = [-5 5];
    else                   % Density shifts, in percent
        limits = [-40 40];
    end
    xlim(limits);
    ylim(limits);
    axis('square');
    
    % Plot data
    xVals = sampledConeParams(:,paramIndsToPlot(i));  % Predicted parameters
    for j = 1:length(testWls)
        yVals = recoveredParamsInd(j,paramIndsToPlot(i));
        plot(xVals,yVals,'o ','MarkerSize',5,'MarkerFaceColor',...
            colorString(i),'MarkerEdgeColor',colorString(i));
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
end