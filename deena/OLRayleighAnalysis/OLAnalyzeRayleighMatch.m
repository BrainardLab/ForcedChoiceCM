function OLAnalyzeRayleighMatch(subjID,sessionNums,varargin)
% Function for analyzing Rayleigh matching data from live subjects.
% Syntax:
%   OLAnalyzeRayleighMatch(subjID, sessionNums)
%
% Description:
%    Analyzes data from human OneLight Rayleigh matching experiments, and
%    uses results to estimate cone individudal difference parameters. Takes
%    in a subject ID and a list of session numbers to analyze. Collects the
%    data from the various sessions, produces a Pitt diagram, and
%    computes
%    the mean and standard deviation of matches for each set of wavelengths.
%    Then, averages the radiometer-measured match spds for each set of
%    wavelengths, and uses these averaged spds to estimate cone individual
%    difference parameters. There is also optional code to recreate matches
%    with a simulated observer as an estimate of observer noise.
%
%    We assume that all sessions share the same adjustment length and that
%    the same number of matches were made for each reference wavelength 
%
% Inputs:
%    subjID              - Character vector of subject ID
%    sessionNums         - Numeric vector of desired session numbers to
%                          analyze.
% Outputs:
%    None (saves a file and figures)
%
% Optional key-value pairs:
%    'LMEqualOD'      -Logical. If true, constrains the L and M cone optical
%                      densities to be equal during cone parameter search.
%                      Default is false.
%    'dlens0'         -Logical. If true, constrains the lens pigment density
%                      to be 0 during cone parameter search. Default is true.
%    'dmac0'          -Logical. If true, constrains the macular pigment
%                      density to be 0 during cone parameter search.
%                      Default is true.
%    'OD0'            -Logical. If true, constrains the optical densities to 
%                      be 0 during cone parameter search. Default is false.
%    'lambdaMax0'     -Logical. If true, constrains lambda max parameters to 
%                      be 0. Default is false.
%    'S0'             -Logical. If true, constrains S lambda max and optical
%                      density (params 5 and 8) to be 0. Default is true.
%    'errScalar'      -Integer for scaling the match error in parameter
%                      fits, to improve search. Default is 100.
%    'restrictBySD'   -Logical. If true, adds lower and upper bounds on all
%                      paramters to keep them within a specified number of 
%                      standard deviations of their means during 
%                      optimization. Default is true.
%    'avgSpds'        -Logical. If true, fits cone parameters based on
%                      averaged spds for each reference wavelength. Default
%                      is false.
%    'multipleParamFits' -Logical. If true, fits parameters using several
%                         different procedures, and makes a cone contrast
%                         plot for each. Note that this overrides 
%                         constraint settings for the lambda max and OD 
%                         params that are entered as key-value pairs. 
%                         Default is true.
%    'makeBarPlots'   -Logical. If true, makes bar plots of cone
%                      excitations for averaged spds. Default is false.
%    'plotPredictedMatches' -Logical. If true, plots predicted and observed
%                            matches for the fit observer. Default is false.
%    'minimizeConeErr'-Logical. If true, minimizes cone exictation error
%                      instead of opponent contrast difference. Default 
%                      is false.
%    'checkOrderEffect'-Logical. If true, the cone excitation plot
%                       highlights which trials were run with the primary
%                       mixture first, and which with the reference first.
%                       Default is false.
%    'sdDensity'     -Number of allowed standard deviations for density 
%                     parameters (1:5) in fit. Default is 3.
%    'sdLambdaMax'   -Number of allowed standard deviations for lambda max
%                     parameters (6:8) in fit. Default is 3.
%    'nCrossValRuns' -Number of times to run the overall cross-validation
%                     procedure. Default is 10.

% History:
%   2/19/21  dce       Wrote it.
%   3/16/21  dce       Added plot saving
%   4/22/21  dce       Added plots for individual excitations
%   4/28/21  dce       Switched to not average spds in fit, restructed and
%                      edited
%   05/09/21  dce      Added option to fit params using cone excitation 
%                      difference instead of opponent contrast.
%   06/04/21  dce      Changed to reflect edits to OLRayleighMatch output
%                      file structure and to findObserverParameters options
%   06/09/21  dce      Added option to check for order effects in cone 
%                      response figures (highlight whether primary or ref 
%                      was shown first)
%   06/16/21  dce      Added an option to do several different cone fits, 
%                      added analysis of deviation from predicted match, 
%                      removed some summary stats 
%   06/21/21  dce      Fixed positioning of reversals, and changed which
%                      spectra are used to calculate match
%   07/05/21  dce      Added option to set number of sds in parameter fit, 
%                      added fits with varying lens density.

close all; 
% Parse input
p = inputParser;
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('OD0',false,@(x)(islogical(x)));
p.addParameter('lambdaMax0',false,@(x)(islogical(x)));
p.addParameter('S0',true,@(x)(islogical(x)));
p.addParameter('sdDensity',3,@(x)(isnumeric(x)));
p.addParameter('sdLambdaMax',3,@(x)(isnumeric(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('avgSpds',false,@(x)(islogical(x)));
p.addParameter('minimizeConeErr',false,@(x)(islogical(x)));
p.addParameter('makeBarPlots',false,@(x)(islogical(x)));
p.addParameter('plotPredictedMatches',false,@(x)(islogical(x)));
p.addParameter('checkOrderEffect',false,@(x)(islogical(x)));
p.addParameter('multipleParamFits',true,@(x)(islogical(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.addParameter('nCrossValRuns',10,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Error checking 
if p.Results.estNoise && ~p.Results.makePittDiagram
    error('Set "make Pitt diagram" to true to estimate noise');
end 

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
refFirst = [];         % Was the reference light shown first for a given match?

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
    
    % Add session data to collected data - light wavelength information
    % and primary ratios/ref intensities
    sessionData = load(outputFile);
    lightCombos = [lightCombos;sessionData.lightCombosFull];
    
    if size(sessionData.primaryRatios)==[2 1]
        primaryRatios = [primaryRatios;sessionData.primaryRatios];
        refIntensities = [refIntensities;sessionData.testIntensities];
    else
        primaryRatios = [primaryRatios;sessionData.primaryRatios'];
        refIntensities = [refIntensities;sessionData.testIntensities'];
    end 
    % Find scale factors for spds
    p1Scales = [p1Scales;sessionData.p1Scale(sessionData.testWls==lightCombos(i,3))];
    p2Scales = [p2Scales;sessionData.p2Scale(sessionData.testWls==lightCombos(i,3))];
    refScales = [refScales;sessionData.testScale(sessionData.testWls==lightCombos(i,3))];
    
    % Collect data for individual matches 
    for j = 1:sessionData.nObserverMatches
        % Find match spds (predicted)
        fName = fullfile(outputDir,[subjID '_' num2str(sessionNums(i))...
            '_' num2str(j) '.mat']);
        trialData = load(fName);
        [rSpds,pSpds,~,~] = getMatchData(fName,'averageSpds',false);
        spdLength = size(rSpds,1);
        predRefSpds = [predRefSpds,reshape(rSpds,[spdLength numel(rSpds)/spdLength])];
        predPrimarySpds = [predPrimarySpds,reshape(pSpds,[spdLength numel(pSpds)/spdLength])];
        trialData.staircaseTestFirst = [true false]; % Debugging - remove
        refFirst = [refFirst,trialData.staircaseTestFirst];
    end
    
    % Add radiometer data to collected data
    radiometerData = load(measFile);
    measPrimarySpds = [measPrimarySpds,radiometerData.measuredPrimarySpds'];
    measRefSpds = [measRefSpds,radiometerData.measuredRefSpds'];
    darkSpds = [darkSpds, repmat(radiometerData.measuredDarkSpd',1,length(sessionData.testIntensities))];
end
% Duplicate light combos array if interleaved
if sessionData.interleaveStaircases
    lightCombos = repelem(lightCombos,2,1);
    p1Scales = repelem(p1Scales,2);
    p2Scales = repelem(p2Scales,2);
    refScales = repelem(refScales,2);
end 

%% Sort collected data by unique wavelength combo, and find average spds
% Identify the different sets of wavelengths used
[matchWls,uniqueWlIndices] = unique(lightCombos,'rows'); % Unique wavelengths tested
[nMatchWls,~] = size(matchWls);                          % Number of unique wls tested
[spdLength,~] = size(measRefSpds);                       % Spd length
nRepeats = size(measRefSpds,2)/nMatchWls;                % Number of times each match was repeated

% Extract data for each set of wavelengths, and average spds.
meanPrimarySpds= zeros(spdLength,nMatchWls);
meanRefSpds = zeros(spdLength,nMatchWls);
measPrimarySpdsByWl = zeros(spdLength,nRepeats,nMatchWls);
measRefSpdsByWl = zeros(spdLength,nRepeats,nMatchWls);
for i = 1:nMatchWls
    measPrimarySpdsByWl(:,:,i) = measPrimarySpds(:,all(lightCombos==matchWls(i,:),2));
    meanPrimarySpds(:,i) = mean(measPrimarySpdsByWl(:,:,i),2);
    
    measRefSpdsByWl(:,:,i) = measRefSpds(:,all(lightCombos==matchWls(i,:),2));
    meanRefSpds(:,i) = mean(measRefSpdsByWl(:,:,i),2);
end

%% Estimate cone fundamentals, and perform cross-validation if we are 
% performing multiple versions of the fit.
% Are we using averaged spds or not?
if p.Results.avgSpds
    primarySpdsFit = meanPrimarySpds;
    refSpdsFit = meanRefSpds;
else
    primarySpdsFit = measPrimarySpds;
    refSpdsFit = measRefSpds;
end

% Generate standard observer
stdObs = genRayleighObserver('age',sessionData.age,'fieldSize',...
    sessionData.fieldSize,'opponentParams',sessionData.opponentParams,...
    'coneVec',zeros(1,8));

% Fit cone fundamentals
if ~p.Results.multipleParamFits
    % Perform one version of the fit, as specified by key-value pairs
    estConeParams = zeros(2,8);
    estObs = {stdObs,[]};
    nConeParams = 1;
    [estConeParams(2,:),~,estObs{2}] = findObserverParameters(refSpdsFit,primarySpdsFit,...
        'age',sessionData.age,'fieldSize',sessionData.fieldSize,...
        'opponentParams',sessionData.opponentParams,'dlens0',p.Results.dlens0,...
        'dmac0',p.Results.dmac0,'LMEqualOD',p.Results.LMEqualOD,...
        'restrictBySd',p.Results.restrictBySd,'OD0',p.Results.OD0,...
        'S0',p.Results.S0,'lambdaMax0',p.Results.lambdaMax0,'minimizeConeErr',...
        p.Results.minimizeConeErr,'sdDensity',p.Results.sdDensity,...
        'sdLambdaMax',p.Results.sdLambdaMax,'errScalar',p.Results.errScalar);
else % Perform several versions of the fit, and cross-validation
    
    % Define bounds. The limits are entered in matrix rows in the following 
    % order: standard, unconstrained, lock OD, lock lambda max, LM equal
    % OD, unconstrained+vary lens,lock OD + vary lens, lock lambda max + 
    % vary lens, LM equal OD + vary lens. 
    nConeParams = 9;
    estObs = cell(1,nConeParams);
    sds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3]; % Parameter standard deviations
    scaleFactors = [repmat(p.Results.sdDensity,1,5), repmat(p.Results.sdLambdaMax,1,3)];
   
    % Manually enter variable parameters and limits for each model
    lowerBounds = zeros(9,8);
    lowerBounds(2,:) = [0 0 1 1 0 1 1 0].*sds.*scaleFactors*-1;
    lowerBounds(3,:) = [0 0 0 0 0 1 1 0].*sds.*scaleFactors*-1;
    lowerBounds(4,:) = [0 0 1 1 0 0 0 0].*sds.*scaleFactors*-1;
    lowerBounds(5,:) = [0 0 1 1 0 1 1 0].*sds.*scaleFactors*-1;   
    lowerBounds(6,:) = [1 0 1 1 0 1 1 0].*sds.*scaleFactors*-1;
    lowerBounds(7,:) = [1 0 0 0 0 1 1 0].*sds.*scaleFactors*-1;
    lowerBounds(8,:) = [1 0 1 1 0 0 0 0].*sds.*scaleFactors*-1;
    lowerBounds(9,:) = [1 0 1 1 0 1 1 0].*sds.*scaleFactors*-1;
    
    upperBounds = zeros(9,8);
    upperBounds(2,:) = [0 0 1 1 0 1 1 0].*sds.*scaleFactors;
    upperBounds(3,:) = [0 0 0 0 0 1 1 0].*sds.*scaleFactors;
    upperBounds(4,:) = [0 0 1 1 0 0 0 0].*sds.*scaleFactors;
    upperBounds(5,:) = [0 0 1 1 0 1 1 0].*sds.*scaleFactors;
    upperBounds(6,:) = [1 0 1 1 0 1 1 0].*sds.*scaleFactors;
    upperBounds(7,:) = [1 0 0 0 0 1 1 0].*sds.*scaleFactors;
    upperBounds(8,:) = [1 0 1 1 0 0 0 0].*sds.*scaleFactors;
    upperBounds(9,:) = [1 0 1 1 0 1 1 0].*sds.*scaleFactors;
    
    % Equality constraint for LM equal OD
    AEq = cell(1,nConeParams);
    BEq = cell(1,nConeParams);
    AEq{5} = [0 0 1 -1 0 0 0 0];
    AEq{9} = [0 0 1 -1 0 0 0 0];
    BEq{5} = 0; 
    BEq{9} = 0; 
    
    % Perform fits
    estConeParams = zeros(nConeParams,8);
    for mm = 1:nConeParams
        [estConeParams(mm,:),~,estObs{mm}] = ...
            findObserverParameters(refSpdsFit,primarySpdsFit,...
            'age',sessionData.age,'fieldSize',sessionData.fieldSize,...
            'opponentParams',sessionData.opponentParams,'initialConeParams',...
            zeros(1,8),'minimizeConeErr',p.Results.minimizeConeErr,...
            'lowerBounds',lowerBounds(mm,:),'upperBounds',upperBounds(mm,:),...
            'AEq',AEq{mm},'BEq',BEq{mm},'errScalar',p.Results.errScalar);
    end
    
    % Run cross-validation program 
    modelCrossValError = ...
    crossValidateRayleighMatch(measPrimarySpdsByWl,measRefSpdsByWl,...
    lowerBounds,upperBounds,AEq,BEq,p.Results.nCrossValRuns,'errScalar',...
    p.Results.errScalar,'age',sessionData.age,'fieldSize',sessionData.fieldSize,...
    'initialConeParams',zeros(1,8),'opponentParams',sessionData.opponentParams);

    % Make bar plot of error  
    crossValErrPlot = figure();
    bar(modelCrossValError);
    plotNames = {'Standard','Unconstrained','Lock OD',...
        'Lock Lambda Max','Equal LM OD','Unconstrained+1','Lock OD+1',...
        'Lock Lambda Max+1','Equal LM OD+1'};
    set(gca,'xticklabel',plotNames);
    text(1:9,modelCrossValError,num2str(modelCrossValError','%0.2f'),...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
    title([subjID ' Cross Validated Fit Error'],'interpreter','none');
    crossValErrPlot.Position = [100 100 1000 400];
    NicePlot.exportFigToPDF([resFile '_crossValErr.pdf'],...
        crossValErrPlot,300); 
end 

%% Make cone response figure for each set of fit params
plotColors = 'rkbgcmrkbgcm';
refWls = unique(lightCombos(:,3));

if p.Results.multipleParamFits
    plotFNames = {'_standard','_unconstrained','_lockOD','_lockLambdaMax',...
        '_LMEqualOD','_unconstrainedLens','_lockODLens','_lockLambdaMaxLens',...
        '_LMEqualODLens'};
else 
    plotNames = {''};
    plotFNames = {''};
end 

% Data arrays to store cone responses
primaryConeRes = cell(nConeParams,nMatchWls);
refConeRes = cell(nConeParams,nMatchWls);
primaryLMinusM = cell(nConeParams,nMatchWls);
refLMinusM = cell(nConeParams,nMatchWls);
primaryLPlusM = cell(nConeParams,nMatchWls);
refLPlusM = cell(nConeParams,nMatchWls);

% Loop through each set of cone parameters, and make a plot
for kk = 1:nConeParams
    legendHandles = [];
    legendEntries = {};
    
    fitConeDiffPlot = figure();
    hold on;
    xlim([0 0.35]);
    ylim([0.015 0.035]);
    xlabel('(L - M)/(L+M)');
    ylabel('L + M');
    title([subjID ' Cone Response Difference ' plotNames{kk}],'interpreter','none');
        
    for i = 1:nMatchWls
        % Extract relevant spds
        sessionInds = all(lightCombos==matchWls(i,:),2);
        measPrimarySpdsTrial = measPrimarySpds(:,sessionInds);
        measRefSpdsTrial = measRefSpds(:,sessionInds);
        
        % Extract relevant refFirst data
        trialRefFirst = refFirst(sessionInds);
        trialRefFirst = logical(trialRefFirst);
        
        % Compute and store excitations
        primaryConeRes{kk,i} = estObs{kk}.T_cones*measPrimarySpdsTrial;
        refConeRes{kk,i} = estObs{kk}.T_cones*measRefSpdsTrial;
        primaryLMinusM{kk,i} = primaryConeRes{kk,i}(1,:)-primaryConeRes{kk,i}(2,:);
        primaryLPlusM{kk,i} = primaryConeRes{kk,i}(1,:)+primaryConeRes{kk,i}(2,:);
        refLMinusM{kk,i} = refConeRes{kk,i}(1,:)-refConeRes{kk,i}(2,:);
        refLPlusM{kk,i} = refConeRes{kk,i}(1,:)+refConeRes{kk,i}(2,:);
     
        % Add points to plot
        figure(fitConeDiffPlot);
        hold on;
        a = plot(primaryLMinusM{kk,i}./primaryLPlusM{kk,i},primaryLPlusM{kk,i},...
            [plotColors(i) '* ']);
        if p.Results.checkOrderEffect
            plot(refLMinusM{kk,i}(trialRefFirst)./refLPlusM{kk,i}(trialRefFirst),...
                refLPlusM{kk,i}(trialRefFirst),[plotColors(i) 'o ']);
            plot(refLMinusM{kk,i}(~trialRefFirst)./refLPlusM{kk,i}(~trialRefFirst),...
                refLPlusM{kk,i}(~trialRefFirst),[plotColors(i) 's ']);
        else
            plot(refLMinusM{kk,i}./refLPlusM{kk,i},refLPlusM{kk,i},[plotColors(i) 'o ']);
        end
        
        % Highlight second session points
        if length(refLMinusM{kk,i}) >=4
            if p.Results.checkOrderEffect
                selectionArr = logical([0 0 refFirst(3:4)]);
                selectionArr2 = logical([0 0 ~refFirst(3:4)]);
                plot(refLMinusM{kk,i}(selectionArr)./refLPlusM{kk,i}(selectionArr),refLPlusM{kk,i}(selectionArr),...
                    'yo','MarkerFaceColor','Yellow','MarkerSize',5);
                plot(refLMinusM{kk,i}(selectionArr2)./refLPlusM{kk,i}(selectionArr2),refLPlusM{kk,i}(selectionArr2),...
                    'ys','MarkerFaceColor','Yellow','MarkerSize',5);
            else
                plot(refLMinusM{kk,i}(3:4)./refLPlusM{kk,i}(3:4),refLPlusM{kk,i}(3:4),...
                    'yo','MarkerFaceColor','Yellow','MarkerSize',5);
            end
        end
        
        % Add lines connecting each pair of primary/reference points
        for j = 1:length(primaryConeRes{kk,i}(1,:))
            plot([primaryLMinusM{kk,i}(j)/primaryLPlusM{kk,i}(j)...
                refLMinusM{kk,i}(j)/refLPlusM{kk,i}(j)],...
                [primaryLPlusM{kk,i}(j) refLPlusM{kk,i}(j)],'k-');
        end
        legendHandles = [legendHandles,a];
        legendEntries{end+1} = num2str(refWls(i));
    end
    
    % Add legends and explanatory text labels, and save figures
    if p.Results.checkOrderEffect
        txtLabel = {'Yellow = second session','Square = primary first'};
    else
        txtLabel = 'Yellow = second session';
    end
    figure(fitConeDiffPlot);
    legend(legendHandles,legendEntries);
    text(0.05,0.03,txtLabel);
    NicePlot.exportFigToPDF([resFile '_fitConeDiffs' plotFNames{kk} '.pdf'],...
        fitConeDiffPlot,300);
end

%% Make plots of deviation from nominal match (optional)
% Loop through each set of fit parameters, and make a separate figure 
primarySpdsPredictedMatch = cell(1,nConeParams);
refSpdsPredictedMatch = cell(1,nConeParams);
primaryConeResPredictedMatch = cell(1,nConeParams);
refConeResPredictedMatch = cell(1,nConeParams);
refLMinusMPredicted = cell(1,nConeParams);
refLPlusMPredicted = cell(1,nConeParams);
primaryLMinusMPredicted = cell(1,nConeParams);
primaryLPlusMPredicted = cell(1,nConeParams);

if p.Results.plotPredictedMatches
    % Preallocate arrays
    for kk = 1:nConeParams
        primarySpdsPredictedMatch{kk} = zeros(spdLength,nMatchWls);
        refSpdsPredictedMatch{kk} = zeros(spdLength,nMatchWls);
        primaryConeResPredictedMatch{kk} = zeros(3,nMatchWls);
        refConeResPredictedMatch{kk} = zeros(3,nMatchWls);
        primaryLMinusMPredicted{kk} = zeros(1,nMatchWls);
        refLMinusMPredicted{kk} = zeros(1,nMatchWls);
        primaryLPlusMPredicted{kk} = zeros(1,nMatchWls);
        refLPlusMPredicted{kk} = zeros(1,nMatchWls);
    end
    % Find predicted matches for each fit observer at each set of wavelengths tested
    for i = 1:nMatchWls
        lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
            sessionData.adjustmentLength,matchWls(i,1),matchWls(i,2),matchWls(i,3),...
            p1Scales(uniqueWlIndices(i)),p2Scales(uniqueWlIndices(i)),...
            refScales(uniqueWlIndices(i)));
        lightData = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
            'precomputedStartStops',lightFile));
        for kk = 1:nConeParams
            [rSpdsPredictedMatch,pSpdsPredictedMatch,~,~] = ...
                searchPredictedRayleighMatch(lightData.testSpdsPredicted,...
                lightData.primarySpdsPredicted,estObs{kk});
            primarySpdsPredictedMatch{kk}(:,i) = pSpdsPredictedMatch;
            refSpdsPredictedMatch{kk}(:,i) = rSpdsPredictedMatch;
            primaryConeResPredictedMatch{kk}(:,i) = estObs{kk}.T_cones*pSpdsPredictedMatch;
            refConeResPredictedMatch{kk}(:,i) = estObs{kk}.T_cones*rSpdsPredictedMatch;
            
            primaryLMinusMPredicted{kk}(i) =...
                primaryConeResPredictedMatch{kk}(1,i)-primaryConeResPredictedMatch{kk}(2,i);
            primaryLPlusMPredicted{kk}(i) = ...
                primaryConeResPredictedMatch{kk}(1,i)+primaryConeResPredictedMatch{kk}(2,i);
            refLMinusMPredicted{kk}(i) = ...
                refConeResPredictedMatch{kk}(1,i)-refConeResPredictedMatch{kk}(2,i);
            refLPlusMPredicted{kk}(i) = ...
                refConeResPredictedMatch{kk}(1,i)+refConeResPredictedMatch{kk}(2,i);
        end
    end
    for kk = 1:nConeParams
        coneDeviationFig = figure();
        hold on;
        legendHandles = [];
        legendEntries = {};
        xlabel('(L-M)/(L+M)');
        ylabel('L+M');
        title([subjID ' Deviation From Predicted Match- ' plotNames{kk}],'interpreter','none');
        txtLabel = {'yellow = predicted matches'};
        text(0.05,0.03,txtLabel);
        
        % Plot predicted matches (yellow)
        plot(primaryLMinusMPredicted{kk}./primaryLPlusMPredicted{kk},primaryLPlusMPredicted{kk},...
            'y* ');
        plot(refLMinusMPredicted{kk}./refLPlusMPredicted{kk},refLPlusMPredicted{kk},...
            'yo ','MarkerFaceColor','Yellow');
        for j = 1:length(primaryConeResPredictedMatch{kk}(1,:))
            plot([primaryLMinusMPredicted{kk}(j)./primaryLPlusMPredicted{kk}(j),...
                refLMinusMPredicted{kk}(j)./refLPlusMPredicted{kk}(j)],...
                [primaryLPlusMPredicted{kk}(j) refLPlusMPredicted{kk}(j)],'k-');
        end
        
        % Plot observed matches
        for i = 1:nMatchWls
            a = plot(primaryLMinusM{kk,i}./primaryLPlusM{kk,i},primaryLPlusM{kk,i},...
                [plotColors(i) '* ']);
            plot(refLMinusM{kk,i}./refLPlusM{kk,i},refLPlusM{kk,i},...
                [plotColors(i) 'o ']);
            for j = 1:length(primaryConeRes{kk,i}(1,:))
                plot([primaryLMinusM{kk,i}(j)/primaryLPlusM{kk,i}(j)...
                    refLMinusM{kk,i}(j)/refLPlusM{kk,i}(j)],...
                    [primaryLPlusM{kk,i}(j) refLPlusM{kk,i}(j)],'k-');
            end
            legendHandles = [legendHandles,a];
            legendEntries{end+1} = num2str(refWls(i));
        end
        legend(legendHandles,legendEntries)
        NicePlot.exportFigToPDF([resFile '_fitConeMatchDeviations' plotFNames{kk} '.pdf'],...
            coneDeviationFig,300);
    end
end

%% Make cone excitations bar graphs for averaged spds (Optional)
% If multiple fits were performed, this figure is made only for the
% unconstrained fit (without lens excitation)
if p.Results.makeBarPlots
    for i = 1:nMatchWls
        primaryRes = estObs{2}.T_cones*meanPrimarySpds(:,i);
        refRes = estObs{2}.T_cones*meanRefSpds(:,i);
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

%% Save data
save([resFile '_analysis.mat'],'p','lightCombos','primaryRatios','refIntensities',...
    'measPrimarySpds','measRefSpds','predPrimarySpds','predRefSpds','darkSpds',...
    'p1Scales','p2Scales','refScales','meanPrimarySpds','meanRefSpds',...
    'estConeParams','estObs','stdObs','radiometerPrimaryRatios','radiometerRefIntensities',...
    'nominalPrimaryRatios','nominalRefIntensities','nominalPrimaryRatiosStd',...
    'nominalRefIntensitiesStd','matchWls','uniqueWlIndices','primaryConeRes','refConeRes',...
    'primaryConeResStd','refConeResStd','primarySpdsPredictedMatch',...
    'refSpdsPredictedMatch','primaryConeResPredictedMatch','refConeResPredictedMatch',...
    'primaryLMinusM','refLMinusM','primaryLPlusM','refLPlusM');
end