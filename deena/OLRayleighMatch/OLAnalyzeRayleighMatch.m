function [estConeParams,meanPrimarySettings,meanRefSettings,...
    stdErrPrimarySettings,stdErrRefSettings]= ...
    OLAnalyzeRayleighMatch(subjID,sessionNums,varargin)
% Function for analyzing Rayleigh matching data from live subjects.
% Syntax:
%   OLAnalyzeRayleighMatch(subjID, sessionNums)
%
% Description:
%    Analyzes data from human OneLight Rayleigh matching experiments, and
%    uses results to estimate cone individudal difference parameters. Takes
%    in a subject ID and a list of session numbers to analyze. Collects the
%    data from the various sessions, produces a Pitt diagram, and computes
%    the mean and standard deviation of matches for each set of wavelengths.
%    Then, averages the radiometer-measured match spds for each set of
%    wavelengths, and uses these averaged spds to estimate cone individual
%    difference parameters. There is also optional code to recreate matches
%    with a simulated observer as an estimate of observer noise.
%
%    We assume that all sessions share the same adjustment length and the
%    same number of matches per trial. 
%
% Inputs:
%    subjID              - Character vector of subject ID
%    sessionNums         - Numeric vector of desired session numbers to
%                          analyze.
% Outputs:
%    coneParams          - Eight-element vector of cone individual difference
%                          parameters estimated based on subject data. See
%                          findObserverParameters for a full description.
%    meanPrimarySettings - Vector containing the average primary mixture
%                          settings, with one entry for each combination of
%                          primary/reference wavelengths tested.
%    meanRefSettings     - Vector containing the average reference
%                          settings, with one entry for each combination of
%                          primary/reference wavelengths tested.
%    stdErrPrimarySettings - Vector containing the standard error of primary
%                          mixture settings, with one entry for each
%                          combination of primary/reference wavelengths.
%    stdErrRefSettings   - Vector containing the standard error of test
%                          settings, with one entry for each combination
%                          of primary/reference wavelengths tested.
%
% Optional key-value pairs:
%    'estNoise'       -Logical. If true, the program runs simulations with
%                      the estimated cone fundamentals at varying noise
%                      levels, computing the standard error of matches for
%                      each set of primary/reference wavelengths. Default
%                      is false.
%    'LMEqualOD'      -Logical. If true, constrains the L and M cone optical
%                      densities to be equal during cone parameter search.
%                      Default is false.
%    'dlens0'         -Logical. If true, constrains the lens pigment density
%                      to be 0 during cone parameter search. Default is true.
%    'dmac0'          -Logical. If true, constrains the macular pigment 
%                      density to be 0 during cone parameter search. 
%                      Default is true.
%    'dLM0'           -Logical. If true, constrains the L and M optical
%                      densities to be 0 during cone parameter search. 
%                      Default is false.
%    'restrictBySD'   -Logical. If true, adds lower and upper bounds on all
%                      paramters to keep them within three standard
%                      deviations of their means during optimization. 
%                      Default is true. 


% History:
%   2/19/21  dce       Wrote it.

% Parse input
p = inputParser;
p.addParameter('estNoise',false,@(x)(islogical(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('dLM0',false,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Define results directory 
resDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),subjID);
if ~exist(resDir,'dir')
    mkdir(resDir);
end 
resFile = fullfile(resDir,[subjID, '_', num2str(sessionNums(1)), '_'...
    num2str(sessionNums(end))]);

% Initialize arrays
lightCombos = [];      % Combinations of primary/reference lights used
primaryRatios = [];    % R/G ratios of subject matches
refIntensities = [];   % Reference intensities of subject matches
measPrimarySpds = [];  % Primary spds identified as matches, as measured by the radiometer
measRefSpds = [];      % Reference spds identified as matches, as measured by the radiometer
predPrimarySpds = [];  % Primary spds identified as matches (predicted)
predRefSpds = [];      % Reference spds identified as matches (predicted)
darkSpds = [];         % Dark spds
p1Scales = [];         % Scale factors for first primary spd
p2Scales = [];         % Scale factors for second primary spd
refScales = [];        % Scale factors for reference light

% Collect match positions and radiometer data from each file
for i = 1:length(sessionNums)
    % Names of results files
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'matchFiles',[subjID '_' num2str(sessionNums(i))]);
    outputFile = fullfile(outputDir,[subjID '_' num2str(sessionNums(i)) '_summary.mat']);
    measFile = fullfile(outputDir,[subjID '_' num2str(sessionNums(i)) '_meas.mat']);
    
    % Check if files exists
    if ~exist(outputFile,'file') || ~exist(outputFile,'file')
        error('Specified results file does not exist');
    end
    
    % Add session data to collected data
    sessionData = load(outputFile);
    
    % Manually-entered parameters (which are now saved by program)
    sessionData.age = 22;
    sessionData.testWls = [590 620];
%     sessionData.fieldSize = 2;
%     sessionData.adjustmentLength = 201;
%     sessionData.opponentParams = [40.3908  205.7353   62.9590    1.0000];
%     sessionData.p1Scale = [1 1 1 1];
%     sessionData.p2Scale = 0.2*[1 1 1 1];
%     sessionData.testScale = 0.2*[1 1 1 1];
    
    lightCombos = [lightCombos;sessionData.lightCombosFull];
    primaryRatios = [primaryRatios;sessionData.primaryRatios'];
    refIntensities = [refIntensities;sessionData.testIntensities'];   
    
    % Go through each match 
    for j = 1:length(sessionData.testIntensities)
        % Fill in settings 
        p1Scales = [p1Scales;sessionData.p1Scale(find(sessionData.testWls==lightCombos(i,3)))];
        p2Scales = [p2Scales;sessionData.p2Scale(find(sessionData.testWls==lightCombos(i,3)))];
        refScales = [refScales;sessionData.testScale(find(sessionData.testWls==lightCombos(i,3)))];
        
        % Find match spds (as measured by radiometer)
        fName = fullfile(outputDir,[subjID '_' num2str(sessionNums(i))...
            '_' num2str(j) '.mat']);
        [rSpds,pSpds,~,~] = getMatchData(fName,'averageSpds',false,'nominal',false);
        predRefSpds = [predRefSpds,rSpds];
        predPrimarySpds = [predPrimarySpds,pSpds];
    end
    
    
    % Add radiometer data to collected data
    radiometerData = load(measFile);
    
    % Manually-entered parameters which are now saved by my match sequence
    % programs
%     cal = OLGetCalibrationStructure('CalibrationType',getpref('ForcedChoiceCM','currentCal'));
%     [~,c] = size(cal.computed.pr650M);
%     radiometerData.darkSpd = OLPrimaryToSpd(cal,zeros(c,1));
    
    measPrimarySpds = [measPrimarySpds,radiometerData.measuredPrimarySpds'];
    measRefSpds = [measRefSpds,radiometerData.measuredTestSpds'];
    darkSpds = [darkSpds, repmat(radiometerData.measuredDarkSpd',1,length(sessionData.testIntensities))];
end

%% Sort collected data by unique wavelength combo, and compute summary
%% statistics
% Identify the different sets of wavelengths used
matchWls = unique(lightCombos,'rows');
[nMatchWls,~] = size(matchWls);
[spdLength,~] = size(measRefSpds);

% Initialize arrays for storing statistics
meanPrimarySettings = zeros(1,nMatchWls);
meanRefSettings = zeros(1,nMatchWls);
stdErrPrimarySettings = zeros(1,nMatchWls);
stdErrRefSettings = zeros(1,nMatchWls);
meanPrimarySpds= zeros(spdLength,nMatchWls);
meanRefSpds = zeros(spdLength,nMatchWls);

% Extract data for each set of wavelengths, compute summary statistics, and
% average spds.
for i = 1:nMatchWls
    primarySettings = primaryRatios(all(lightCombos==matchWls(i,:),2));
    meanPrimarySettings(i) = mean(primarySettings);
    stdErrPrimarySettings(i) = std(primarySettings)/sqrt(length(primarySettings));
    
    refSettings = refIntensities(all(lightCombos==matchWls(i,:),2));
    meanRefSettings(i) = mean(refSettings);
    stdErrRefSettings(i) = std(refSettings)/sqrt(length(refSettings));
    
    primarySpds = measPrimarySpds(:,all(lightCombos==matchWls(i,:),2));
    meanPrimarySpds(:,i) = mean(primarySpds,2);
    
    refSpds = measRefSpds(:,all(lightCombos==matchWls(i,:),2));
    meanRefSpds(:,i) = mean(refSpds,2);
end

%% Estimate cone fundamentals based on averaged spds
% Use variables from last results file
[estConeParams,~,~] = findObserverParameters(meanRefSpds,meanPrimarySpds,...
    'age',sessionData.age,'fieldSize',sessionData.fieldSize,...
    'opponentParams',sessionData.opponentParams,'dlens0',p.Results.dlens0,...
    'dmac0',p.Results.dmac0,'LMEqualOD',p.Results.LMEqualOD,...
    'restrictBySd',p.Results.restrictBySd,'dLM0',p.Results.dLM0);
estObs = genRayleighObserver('age',sessionData.age,'fieldSize',...
    sessionData.fieldSize,'opponentParams',sessionData.opponentParams,...
    'coneVec',estConeParams);
stdObs = genRayleighObserver('age',sessionData.age,'fieldSize',...
    sessionData.fieldSize,'opponentParams',sessionData.opponentParams,...
    'coneVec',zeros(1,8));

%% Make a cone excitations bar graph
for i = 1:nMatchWls
    primaryRes = estObs.T_cones*meanPrimarySpds(:,i);
    refRes = estObs.T_cones*meanRefSpds(:,i);
    primaryResStd = stdObs.T_cones*meanPrimarySpds(:,i);
    refResStd = stdObs.T_cones*meanRefSpds(:,i);
    
    figure();
    hold on;
    subplot(2,1,1);   % L cone excitations
    bar([primaryRes(1),refRes(1); primaryResStd(1),refResStd(1)]);
    title('L Cone Excitations');
    legend('Primary Light','Reference Light');
    ylim([0 0.03]);
    names ={'Fitted Observer','Standard Observer'};
    set(gca,'xticklabel', names)
    ylabel('Relative Response Intensity');
  
    subplot(2,1,2);   % M cone excitations  
    bar([primaryRes(2),refRes(2);primaryResStd(2),refResStd(2)]);
    title('M Cone Excitations');
    legend('Primary Light','Reference Light');
    ylim([0 0.03]);
    
    names ={'Fitted Observer','Standard Observer'};
    set(gca,'xticklabel', names)
    ylabel('Relative Response Intensity');
    
    theTitle = sprintf('Cone Excitations for Fitted Observers, ref = %g',matchWls(i,3));
    sgtitle(theTitle);
end 

%% Make summary Pitt diagram
% Bring points to a common scaling
commonP1ScaleFactors = p1Scales/max(p1Scales);
commonP2ScaleFactors = p2Scales/max(p2Scales);
if any(commonP1ScaleFactors~=commonP2ScaleFactors)
    error('P1 and P2 must be scaled proportionally')
end
commonRefScaleFactors = refScales/max(refScales);

% Compute Pitt diagram points for radiometer measurements
radiometerPrimaryRatios = [];
radiometerRefIntensities = [];
nominalPrimaryRatios = [];
nominalRefIntensities = [];
nominalPrimaryRatiosStd = [];
nominalRefIntensitiesStd = [];

for i = 1:length(primaryRatios)
    % Define output light file we're using
    lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
        sessionData.adjustmentLength,lightCombos(i,1),lightCombos(i,2),...
        lightCombos(i,3),max(p1Scales),max(p2Scales),max(refScales));
    lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'precomputedStartStops',lightFile);
    lightSettings = load(lightFileName);
    
    % Find Pitt points
    [refIntensityR,pRatioR] = OLSpdToPittPoint(measRefSpds(:,i),...
        measPrimarySpds(:,i),darkSpds(:,i),lightFileName);
    radiometerPrimaryRatios = [radiometerPrimaryRatios,pRatioR];
    radiometerRefIntensities = [radiometerRefIntensities,refIntensityR];
    
    % Find nominal match for standard observer and for estimated observer
    [~,~,rIndex,pIndex] = searchPredictedRayleighMatch(...
        lightSettings.testSpdsPredicted,lightSettings.primarySpdsPredicted,...
        estObs);
    [~,~,rIndexStd,pIndexStd] = searchPredictedRayleighMatch(...
        lightSettings.testSpdsPredicted,lightSettings.primarySpdsPredicted,...
        stdObs);
    
    nominalPrimaryRatios = [nominalPrimaryRatios,lightSettings.p1Scales(pIndex)];
    nominalRefIntensities = [nominalRefIntensities,lightSettings.testScales(rIndex)];   
    nominalPrimaryRatiosStd = [nominalPrimaryRatiosStd,lightSettings.p1Scales(pIndexStd)];
    nominalRefIntensitiesStd = [nominalRefIntensitiesStd,lightSettings.testScales(rIndexStd)];
end

pittPlot = figure();
hold on;
xlabel('Primary Ratio');
ylabel('Test Intensity');
plot(primaryRatios.*commonP1ScaleFactors,refIntensities.*commonRefScaleFactors,...
    'bs ','MarkerFaceColor','b','MarkerSize',5);
plot(radiometerPrimaryRatios,radiometerRefIntensities,...
    'r* ','MarkerSize',5);
plot(nominalPrimaryRatios,nominalRefIntensities,...
    'go ','MarkerSize',5,'MarkerFaceColor','g');
plot(nominalPrimaryRatiosStd,nominalRefIntensitiesStd,...
    'mo ','MarkerSize',5,'MarkerFaceColor','m');
legend('Settings','Radiometer Spd Fits','Nominal Settings (Estimated Observer)',...
    'Nominal Settings (Standard Observer)');
title([subjID ' Rayleigh Matches'],'interpreter','none');
xlim([0 1]);
ylim([0 1]);

% Save figure 
NicePlot.exportFigToPDF([resFile '_pittPlot.pdf'],...
    pittPlot,300);

%% If desired, conduct matches with a simulated observer at different noise levels
% We simulate using both recovered cone parameters and standard cone
% parameters
% Results are plotted in a separate Pitt diagram 

% Define output arrays 
    
% Find the number of times each match was repeated.
% We assume an equal number of trials were done for each match, and an
% equal number of trials were done in each session.
nRepeats = length(sessionNums)*length(sessionData.primaryRatios)/nMatchWls;

refIntensitiesSim = zeros(nRepeats*nMatchWls,4);
primaryRatiosSim = zeros(nRepeats*nMatchWls,4);
refIntensitiesSimStd = zeros(nRepeats*nMatchWls,4);
primaryRatiosSimStd = zeros(nRepeats*nMatchWls,4);

if p.Results.estNoise
    % Load a Rayleigh match file from the real observer to extract needed
    % experimental parameters
    sampleSessionFile = fullfile(outputDir, [subjID '_'...
        num2str(sessionNums(end)) '_1.mat']);
    sampleSession = load(sampleSessionFile);
    
    % Set preference for simulation to dev file so that we can write from
    % another device
    baseDir = '/Users/deena/Dropbox (Aguirre-Brainard Lab)';
    setpref('ForcedChoiceCM','rayleighDataDir',...
        fullfile(baseDir,'MELA_datadev','Experiments','ForcedChoiceCM','OLRayleighMatch'));
    
    % Run a simulated experiment with the model observer
    for noiseSd = 1:4
        subjIDSim = [subjID '_sim_' num2str(noiseSd)];
        subjIDSimStd = [subjID '_simStd_' num2str(noiseSd)];
        method = 'forcedChoice';
        
        % Note that getMatchSeries computes matches for all possible
        % combinations of the input wavelengths, not necessarily the specific
        % ones that were done here. This should not be relevant if we are only
        % varying reference wavelength.
        [~,~,refIntensitiesSim(:,noiseSd),primaryRatiosSim(:,noiseSd)] = ...
            getMatchSeries(subjIDSim,estConeParams,sessionData.opponentParams,...
            unique(matchWls(:,1))',unique(matchWls(:,2))',unique(matchWls(:,3))',...
            method,'fieldSize',sessionData.fieldSize,'age',sessionData.age,...
            'nObserverMatches',nRepeats,'rayleighPlots',false,'noiseScaleFactor',...
            noiseSd,'averageSpds',false,'p1Scale',sampleSession.p1Scale*[1 1],...
            'p2Scale',sampleSession.p2Scale*[1 1],'testScale',...
            sampleSession.testScale*[1 1],'nReversals',...
            sampleSession.nReversals,'adjustmentLength',...
            sampleSession.adjustmentLength);
        [~,~,refIntensitiesSimStd(:,noiseSd),primaryRatiosSimStd(:,noiseSd)] = ...
            getMatchSeries(subjIDSimStd,zeros(1,8),sessionData.opponentParams,...
            unique(matchWls(:,1))',unique(matchWls(:,2))',unique(matchWls(:,3))',...
            method,'fieldSize',sessionData.fieldSize,'age',sessionData.age,...
            'nObserverMatches',nRepeats,'rayleighPlots',false,'noiseScaleFactor',...
            noiseSd,'averageSpds',false,'p1Scale',sampleSession.p1Scale*[1 1],...
            'p2Scale',sampleSession.p2Scale*[1 1],'testScale',...
            sampleSession.testScale*[1 1],'nReversals',...
            sampleSession.nReversals,'adjustmentLength',sampleSession.adjustmentLength);
    end
    % Reset preference to MELA_data
    baseDir = '/Users/deena/Dropbox (Aguirre-Brainard Lab)';
    setpref('ForcedChoiceCM','rayleighDataDir',...
        fullfile(baseDir,'MELA_data','Experiments','ForcedChoiceCM','OLRayleighMatch'));
    
    refIntensitiesSim = reshape(refIntensitiesSim,nRepeats,nMatchWls,4);
    primaryRatiosSim = reshape(primaryRatiosSim,nRepeats,nMatchWls,4);
    refIntensitiesSimStd = reshape(refIntensitiesSimStd,nRepeats,nMatchWls,4);
    primaryRatiosSimStd = reshape(primaryRatiosSimStd,nRepeats,nMatchWls,4);
    
    % Make a Pitt plot at an intermediate noise level 
    plotNoiseSd = 3;
    simPittPlot = figure();
    hold on;
    xlabel('Primary Ratio');
    ylabel('Test Intensity');
    a = plot(primaryRatios.*commonP1ScaleFactors,refIntensities.*commonRefScaleFactors,...
        'bs ','MarkerFaceColor','b','MarkerSize',5);
    b = plot(radiometerPrimaryRatios,radiometerRefIntensities,...
        'r* ','MarkerSize',4);
    c = plot(primaryRatiosSim(:,:,plotNoiseSd),refIntensitiesSim(:,:,plotNoiseSd),...
        'go ','MarkerFaceColor','g','MarkerSize',3);
    d = plot(primaryRatiosSimStd(:,:,plotNoiseSd),refIntensitiesSimStd(:,:,plotNoiseSd),...
        'mo ','MarkerFaceColor','m','MarkerSize',2);
    legend([a(1),b(1),c(1),d(1)],'Settings','Radiometer Spd Fits',...
        'Simulated Observer Matches (Estimated Cones)',...
        'Simulated Observer Matches (Standard Cones)');
    title([subjID ' Rayleigh Matches - Simulated Observer Comparison']);
    xlim([0 1]);
    ylim([0 1]);
    
    % Save figure 
    NicePlot.exportFigToPDF([resFile '_simPittPlot.pdf'],...
        simPittPlot,300);
    
    % Compute standard error for matches with each set of primary/reference
    % wavelength
    refIntensitiesSimStdErr = zeros(4,nMatchWls);
    primaryRatiosSimStdErr = zeros(4,nMatchWls);
    for i = 1:4
        for j = 1:nMatchWls
            refIntensitiesSimStdErr(i,j) = std(refIntensitiesSim(:,j,i))/sqrt(length(refIntensitiesSim(:,j,i)));
            primaryRatiosSimStdErr(i,j) = std(primaryRatiosSim(:,j,i))/sqrt(length(primaryRatiosSim(:,j,i)));
        end
    end
end

%% Plot of measured and predicted spds
% Can turn on if desired
plotSpdsMeasPred = false;
if plotSpdsMeasPred
    spdPlot = figure();
    wls = SToWls([380 2 201]);
    [nMatches,~] = size(predRefSpds); 
    for i = 1:nMatches
        subplot(2,4,i);
        hold on;
        b1 = plot(wls,predRefSpds(:,i),'b',wls,predPrimarySpds(:,i),'b','LineWidth',3);
        b2 = plot(wls,measRefSpds(:,i),'r',wls,measPrimarySpds(:,i),'r','LineWidth',1);
        legend([b1(1) b2(1)],'Predicted Spds','Measured Spds');
        xlabel('Wavelength (nm)');
        ylabel('Power');
        theTitle = sprintf('Match %g',i);
        title(theTitle);
    end
    NicePlot.exportFigToPDF([resFile '_spdPlot.pdf'],...
        spdPlot,300);
end
% Save data
save([resFile '_analysis.mat'],'p','lightCombos','primaryRatios',...
    'refIntensities','measPrimarySpds','measRefSpds','predPrimarySpds',...
    'predRefSpds','darkSpds','p1Scales','p2Scales','refScales',...
    'meanPrimarySettings','meanRefSettings','stdErrPrimarySettings',...
    'stdErrRefSettings','meanPrimarySpds','meanRefSpds','commonP1ScaleFactors',...
    'commonP2ScaleFactors','commonRefScaleFactors','estConeParams',...
    'refIntensitiesSim','primaryRatiosSim','refIntensitiesSimStd',...
    'primaryRatiosSimStd');
end