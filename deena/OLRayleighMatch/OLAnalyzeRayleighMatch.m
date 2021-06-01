function OLAnalyzeRayleighMatch(subjID,sessionNums,varargin)
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
%    None (saves a file and figures)
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
%    'avgSpds'        -Logical. If true, fits cone parameters based on
%                      averaged spds for each reference wavelength. Default
%                      is false.
%    'makeBarPlots'   -Logical. If true, makes bar plots of cone
%                      excitations for averaged spds. Default is false.
%    'makePittDiagram'-Logical. If true, makes generalized Pitt diagram.
%                      Default is false.
%    minimizeConeErr  -Logical. If true, minimizes cone exictation error
%                      instead of opponent contrast difference. Default 
%                      is false.


% History:
%   2/19/21  dce       Wrote it.
%   3/16/21  dce       Added plot saving
%   4/22/21  dce       Added plots for individual excitations
%   4/28/21  dce       Switched to not average spds in fit, restructed and
%                      edited
%   05/09/21  dce      Added option to fit params using cone excitation 
%                      difference instead of opponent contrast.
close all; 

% Parse input
p = inputParser;
p.addParameter('estNoise',false,@(x)(islogical(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('dLM0',false,@(x)(islogical(x)));
p.addParameter('dS0',true,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('avgSpds',false,@(x)(islogical(x)));
p.addParameter('minimizeConeErr',false,@(x)(islogical(x)));
p.addParameter('makeBarPlots',false,@(x)(islogical(x)));
p.addParameter('makePittDiagram',false,@(x)(islogical(x)));
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
        'matchFiles',subjID,[subjID '_' num2str(sessionNums(i))]);
    outputFile = fullfile(outputDir,[subjID '_' num2str(sessionNums(i)) '_summary.mat']);
    measFile = fullfile(outputDir,[subjID '_' num2str(sessionNums(i)) '_meas.mat']);
    
    % Check if files exists
    if ~exist(outputFile,'file') || ~exist(outputFile,'file')
        error('Specified results file does not exist');
    end
    
    % Add session data to collected data
    sessionData = load(outputFile);
    
    lightCombos = [lightCombos;sessionData.lightCombosFull];
    primaryRatios = [primaryRatios;sessionData.primaryRatios'];
    refIntensities = [refIntensities;sessionData.testIntensities'];
    
    % Go through each match
    for j = 1:length(sessionData.testIntensities)
        % Fill in scale factors for spds
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
    measPrimarySpds = [measPrimarySpds,radiometerData.measuredPrimarySpds'];
    measRefSpds = [measRefSpds,radiometerData.measuredTestSpds'];
    darkSpds = [darkSpds, repmat(radiometerData.measuredDarkSpd',1,length(sessionData.testIntensities))];
end

%% Sort collected data by unique wavelength combo, and find average spds
% Identify the different sets of wavelengths used
matchWls = unique(lightCombos,'rows');
[nMatchWls,~] = size(matchWls);
[spdLength,~] = size(measRefSpds);

% Extract data for each set of wavelengths and average spds.
meanPrimarySpds= zeros(spdLength,nMatchWls);
meanRefSpds = zeros(spdLength,nMatchWls);
for i = 1:nMatchWls
    primarySpds = measPrimarySpds(:,all(lightCombos==matchWls(i,:),2));
    meanPrimarySpds(:,i) = mean(primarySpds,2);
    
    refSpds = measRefSpds(:,all(lightCombos==matchWls(i,:),2));
    meanRefSpds(:,i) = mean(refSpds,2);
end

%% Estimate cone fundamentals
% Are we using averaged spds or not?
if p.Results.avgSpds
    primarySpdsFit = meanPrimarySpds;
    refSpdsFit = meanRefSpds;
else
    primarySpdsFit = measPrimarySpds;
    refSpdsFit = measRefSpds;
end

% Fit cone fundamentals
[estConeParams,~,~] = findObserverParameters(refSpdsFit,primarySpdsFit,...
    'age',sessionData.age,'fieldSize',sessionData.fieldSize,...
    'opponentParams',sessionData.opponentParams,'dlens0',p.Results.dlens0,...
    'dmac0',p.Results.dmac0,'LMEqualOD',p.Results.LMEqualOD,...
    'restrictBySd',p.Results.restrictBySd,'dLM0',p.Results.dLM0,...
    'dS0',true,'minimizeConeErr',p.Results.minimizeConeErr);

% Standard and fit observers
estObs = genRayleighObserver('age',sessionData.age,'fieldSize',...
    sessionData.fieldSize,'opponentParams',sessionData.opponentParams,...
    'coneVec',estConeParams);
stdObs = genRayleighObserver('age',sessionData.age,'fieldSize',...
    sessionData.fieldSize,'opponentParams',sessionData.opponentParams,...
    'coneVec',zeros(1,8));

%% Make cone excitations bar graphs for averaged spds (Optional)
if p.Results.makeBarPlots
    for i = 1:nMatchWls
        primaryRes = estObs.T_cones*meanPrimarySpds(:,i);
        refRes = estObs.T_cones*meanRefSpds(:,i);
        primaryResStd = stdObs.T_cones*meanPrimarySpds(:,i);
        refResStd = stdObs.T_cones*meanRefSpds(:,i);
        
        avgBarPlot = figure();
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
        
        theTitle = sprintf('Cone Excitations for Fitted Observers, ref = %g (Averaged)',matchWls(i,3));
        sgtitle(theTitle);
        NicePlot.exportFigToPDF([resFile '_' num2str(matchWls(i,3)) '_avgExcitations.pdf'],...
            avgBarPlot,300);
    end
end

%% Make cone response figures
% Setup
plotColors = 'rkbgcmrkbgcm';
refWls = unique(lightCombos(:,3));
legendHandles1 = [];
legendHandles2 = [];
legendEntries = {};

stdConeDiffPlot = figure();
hold on;
xlim([0 0.35]);
ylim([0.015 0.035]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title([subjID ' Cone Response Difference - Standard Cones'],'interpreter','none');

fitConeDiffPlot = figure();
hold on;
xlim([0 0.35]);
ylim([0.015 0.035]);
xlabel('(L - M)/(L+M)');
ylabel('L + M');
title([subjID ' Cone Response Difference - Fit Cones'],'interpreter','none');

% Data arrays
meanPrimaryExcitationsFitL = zeros(1,nMatchWls);
meanRefExcitationsFitL = zeros(1,nMatchWls);
semPrimaryExcitationsFitL = zeros(1,nMatchWls);
semRefExcitationsFitL = zeros(1,nMatchWls);

meanPrimaryExcitationsFitM = zeros(1,nMatchWls);
meanRefExcitationsFitM = zeros(1,nMatchWls);
semPrimaryExcitationsFitM = zeros(1,nMatchWls);
semRefExcitationsFitM = zeros(1,nMatchWls);

meanPrimaryExcitationsStdL = zeros(1,nMatchWls);
meanRefExcitationsStdL = zeros(1,nMatchWls);
semPrimaryExcitationsStdL = zeros(1,nMatchWls);
semRefExcitationsStdL = zeros(1,nMatchWls);

meanPrimaryExcitationsStdM = zeros(1,nMatchWls);
meanRefExcitationsStdM = zeros(1,nMatchWls);
semPrimaryExcitationsStdM = zeros(1,nMatchWls);
semRefExcitationsStdM = zeros(1,nMatchWls);

for i = 1:nMatchWls
    % Extract relevant spds
    sessionInds = all(lightCombos==matchWls(i,:),2);
    measPrimarySpdsTrial = measPrimarySpds(:,sessionInds);
    measRefSpdsTrial = measRefSpds(:,sessionInds);
    
    % Compute excitations
    primaryRes = estObs.T_cones*measPrimarySpdsTrial;
    refRes = estObs.T_cones*measRefSpdsTrial;
    primaryLMinusM = primaryRes(1,:)-primaryRes(2,:);
    primaryLPlusM = primaryRes(1,:)+primaryRes(2,:);
    refLMinusM = refRes(1,:)-refRes(2,:);
    refLPlusM = refRes(1,:)+refRes(2,:);
    
    primaryResStd = stdObs.T_cones*measPrimarySpdsTrial;
    refResStd = stdObs.T_cones*measRefSpdsTrial;
    primaryLMinusMStd = primaryResStd(1,:)-primaryResStd(2,:);
    primaryLPlusMStd = primaryResStd(1,:)+primaryResStd(2,:);
    refLMinusMStd = refResStd(1,:)-refResStd(2,:);
    refLPlusMStd = refResStd(1,:)+refResStd(2,:);
    
    % Compute summary statistics
    meanPrimaryExcitationsFitL(i) = mean(primaryRes(1,:));
    meanRefExcitationsFitL(i) = mean(refRes(1));
    semPrimaryExcitationsFitL(i) = std(primaryRes(1,:))/sqrt(length(primaryRes(1,:)));
    semRefExcitationsFitL(i) = std(refRes(1,:))/sqrt(length(refRes(1,:)));
    
    meanPrimaryExcitationsFitM(i) = mean(primaryRes(2,:));
    meanRefExcitationsFitM(i) =  mean(refRes(2,:));
    semPrimaryExcitationsFitM(i) = std(primaryRes(2,:))/sqrt(length(primaryRes(2,:)));
    semRefExcitationsFitM(i) = std(refRes(2,:))/sqrt(length(refRes(2,:)));
    
    meanPrimaryExcitationsStdL(i) =  mean(primaryResStd(1,:));
    meanRefExcitationsStdL(i) = mean(refResStd(1,:));
    semPrimaryExcitationsStdL(i) = std(primaryResStd(1,:))/sqrt(length(primaryResStd(1,:)));
    semRefExcitationsStdL(i) = std(refResStd(1,:))/sqrt(length(refResStd(1,:)));
    
    meanPrimaryExcitationsStdM(i) = mean(primaryResStd(2,:));
    meanRefExcitationsStdM(i) = mean(refResStd(2,:));
    semPrimaryExcitationsStdM(i) = std(primaryResStd(2,:))/sqrt(length(primaryResStd(2,:)));
    semRefExcitationsStdM(i) = std(refResStd(2,:))/sqrt(length(refResStd(2,:)));
    
    % Add points to plot
    figure(stdConeDiffPlot);
    hold on;
    a  = plot(primaryLMinusMStd./primaryLPlusMStd,...
        primaryLPlusMStd,[plotColors(i) '* ']);
    plot(refLMinusMStd./refLPlusMStd,refLPlusMStd,[plotColors(i) 'o ']);
    if length(refLMinusMStd) >=4
        plot(refLMinusMStd(3:4)./refLPlusMStd(3:4),refLPlusMStd(3:4),...
            'yo','MarkerFaceColor','Yellow','MarkerSize',5);
    end
    legendHandles1 = [legendHandles1,a];
    
    figure(fitConeDiffPlot);
    hold on;
    a = plot(primaryLMinusM./primaryLPlusM,primaryLPlusM,...
        [plotColors(i) '* ']);
    plot(refLMinusM./refLPlusM,refLPlusM,[plotColors(i) 'o ']);
    if length(refLMinusM)>=4
        plot(refLMinusM(3:4)./refLPlusM(3:4),refLPlusM(3:4),...
            'yo','MarkerFaceColor','Yellow','MarkerSize',5);
    end
    legendHandles2 = [legendHandles2,a];
    
    legendEntries{end+1} = num2str(refWls(i));
    
    % Add lines connecting each pair of primary/reference points
    for j = 1:length(primaryRes(1,:))
        figure(stdConeDiffPlot);
        plot([primaryLMinusMStd(j)/primaryLPlusMStd(j)...
            refLMinusMStd(j)/refLPlusMStd(j)],...
            [primaryLPlusMStd(j) refLPlusMStd(j)],'k-');
        
        figure(fitConeDiffPlot);
        plot([primaryLMinusM(j)/primaryLPlusM(j)...
            refLMinusM(j)/refLPlusM(j)],...
            [primaryLPlusM(j) refLPlusM(j)],'k-');
    end
end

% Add legends and save figures
figure(stdConeDiffPlot);
legend(legendHandles1,legendEntries);
NicePlot.exportFigToPDF([resFile '_stdConeDiffs.pdf'],...
    stdConeDiffPlot,300);

figure(fitConeDiffPlot);
legend(legendHandles2,legendEntries);
NicePlot.exportFigToPDF([resFile '_fitConeDiffs.pdf'],...
    fitConeDiffPlot,300);

%% Make summary Pitt diagram
% Compute Pitt diagram points for radiometer measurements, using the
% most extreme scale factors
radiometerPrimaryRatios = [];
radiometerRefIntensities = [];
nominalPrimaryRatios = [];
nominalRefIntensities = [];
nominalPrimaryRatiosStd = [];
nominalRefIntensitiesStd = [];
if p.Results.makePittDiagram
    baseDir = '/Users/deena/Dropbox (Aguirre-Brainard Lab)';
    refWls = unique(lightCombos(:,3));
    for i = 1:length(primaryRatios)
        
        % Define output light file we're using. Either load the data or
        % create file if it does not exist
        lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
            sessionData.adjustmentLength,lightCombos(i,1),lightCombos(i,2),...
            lightCombos(i,3),max(p1Scales),min(p2Scales),max(refScales));
        lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
            'precomputedStartStops',lightFile);
        if ~exist(lightFileName,'file')
            setpref('ForcedChoiceCM','rayleighDataDir',...
                fullfile(baseDir,'MELA_datadev','Experiments','ForcedChoiceCM','OLRayleighMatch'));
            lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
                'precomputedStartStops',lightFile);
            if ~exist(lightFileName,'file')
                OLRayleighMatchLightSettings(lightCombos(i,1),lightCombos(i,2),...
                    lightCombos(i,3),'p1ScaleFactor',max(p1Scales),...
                    'p2ScaleFactor',min(p2Scales),'testScaleFactor',max(refScales),...
                    'adjustmentLength',sessionData.adjustmentLength);
            end
            lightSettings = load(lightFileName);
            setpref('ForcedChoiceCM','rayleighDataDir',...
                fullfile(baseDir,'MELA_data','Experiments','ForcedChoiceCM','OLRayleighMatch'));
        else
            lightSettings = load(lightFileName);
        end
        
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
    ylabel('Reference Intensity');
    legendHandles = [];
    legendEntries = {};
    for i = 1:length(refWls)
        indsForWl = (lightCombos(:,3)==refWls(i));
        pointColor = plotColors(i);
        a = plot(radiometerPrimaryRatios(indsForWl),radiometerRefIntensities(indsForWl),...
            '* ','MarkerSize',5,'MarkerEdgeColor',pointColor);
        plot(nominalPrimaryRatios(indsForWl),nominalRefIntensities(indsForWl),...
            [pointColor 'o '],'MarkerSize',5);
        plot(nominalPrimaryRatiosStd(indsForWl),nominalRefIntensitiesStd(indsForWl),...
            [pointColor 's '],'MarkerSize',5);
        legendHandles = [legendHandles,a];
        legendEntries{end+1} = num2str(refWls(i));
    end
    legend(legendHandles,legendEntries);
    title([subjID ' Rayleigh Matches'],'interpreter','none');
    xlim([0 1]);
    ylim([0 1]);
    
    % Save figure
    NicePlot.exportFigToPDF([resFile '_pittPlot.pdf'],...
        pittPlot,300);
end

%% If desired, conduct matches with a simulated observer at different noise levels
% We simulate using both recovered cone parameters and standard cone
% parameters
% Results are plotted in a separate Pitt diagram

% Find the number of times each match was repeated.
% We assume an equal number of trials were done for each match, and an
% equal number of trials were done in each session.
if p.Results.estNoise
    nRepeats = length(sessionNums)*length(sessionData.primaryRatios)/nMatchWls;
    
    % Define output arrays
    refIntensitiesSim = zeros(nRepeats*nMatchWls,4);
    primaryRatiosSim = zeros(nRepeats*nMatchWls,4);
    refIntensitiesSimStd = zeros(nRepeats*nMatchWls,4);
    primaryRatiosSimStd = zeros(nRepeats*nMatchWls,4);
    
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
    title([subjID ' Rayleigh Matches - Simulated Observer Comparison'],'interpreter','none');
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

%% Plot of measured and predicted spds - not made by default
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

%% Save data
save([resFile '_analysis.mat'],'p','lightCombos','primaryRatios','refIntensities',...
    'measPrimarySpds','measRefSpds','predPrimarySpds','predRefSpds','darkSpds',...
    'p1Scales','p2Scales','refScales','meanPrimarySpds','meanRefSpds',...
    'estConeParams','estObs','stdObs','meanPrimaryExcitationsFitL',...
    'meanPrimaryExcitationsFitM','meanRefExcitationsFitL','meanRefExcitationsFitM',...
    'semPrimaryExcitationsFitL','semPrimaryExcitationsFitM','semRefExcitationsFitL',...
    'semRefExcitationsFitM','meanPrimaryExcitationsStdL',...
    'meanPrimaryExcitationsStdM','meanRefExcitationsStdL','meanRefExcitationsStdM',...
    'semPrimaryExcitationsStdL','semPrimaryExcitationsStdM','semRefExcitationsStdL',...
    'semRefExcitationsStdM','radiometerPrimaryRatios','radiometerRefIntensities',...
    'nominalPrimaryRatios','nominalRefIntensities','nominalPrimaryRatiosStd',...
    'nominalRefIntensitiesStd');
end
