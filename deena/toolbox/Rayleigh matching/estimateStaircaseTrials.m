% This scripts estimates the number of trials used in a staircase Rayleigh
% matching procedure (using getMatchSeries)

% History 
%   11/2/20  dce   Wrote it.

% Rayleigh match settings (can be modified).
nObservers = 20;
coneParamsToVary = [0 0 1 1 0 1 1 0];
baseConeParams = zeros(1,length(coneParamsToVary));
opponentParams = [0.8078 4.1146 1.2592 0.02];
age = 32;
fieldSize = 2;
rayleighPlots = false;
saveResults = false;
adjustmentLength = 3201;
noiseScaleFactor = 0;  % Start with the noiseless case
monochromatic = true;
nObserverMatches = 1;
averageSpds = true;
nReversals = [1 4];
nBelowThreshold = 1;
thresholdScaleFactor = 0.5;
baseSubjID = 'TestNTrials4';

% Observer params
sampledConeParams = sampleRayleighObservers(nObservers,baseConeParams,...
        coneParamsToVary);

% Wavelength params 
p1Wl = 670;
p2Wl = 560;
testWls = [570 590 610 630 650];
p1Scale = 1;
p2Scale = 0.02;
testScale = 0.5;
S = [380 2 201]; 

% Data arrays 
nAdjustmentsFC = zeros(nObservers,1);  
nAdjustmentsThresh = zeros(nObservers,1);

% Loop through observers and compute a series of matches
for i = 1:nObservers
    subjIDFC = [baseSubjID 'FC_' num2str(i)];
    subjIDThreshold = [baseSubjID 'Thresh_' num2str(i)];
    % Compute a series of matches with each method
    getMatchSeries(subjIDFC,sampledConeParams(i,:),opponentParams,...
        p1Wl,p2Wl,testWls,'forcedChoice','age',age,'fieldSize',fieldSize,...
        'monochromatic',monochromatic,'p1Scale',p1Scale,'p2Scale',p2Scale,...
        'testScale',testScale,'sPredicted',S,'rayleighPlots',rayleighPlots,...
        'saveResults',saveResults,'adjustmentLength',adjustmentLength,...
        'noiseScaleFactor',noiseScaleFactor,'averageSpds',averageSpds,...
        'nReversals',nReversals,'nBelowThreshold',nBelowThreshold,...
        'nObserverMatches',nObserverMatches,'thresholdScaleFactor',...
        thresholdScaleFactor);
    getMatchSeries(subjIDThreshold,sampledConeParams(i,:),opponentParams,...
        p1Wl,p2Wl,testWls,'threshold','age',age,'fieldSize',fieldSize,...
        'monochromatic',monochromatic,'p1Scale',p1Scale,'p2Scale',p2Scale,...
        'testScale',testScale,'sPredicted',S,'rayleighPlots',rayleighPlots,...
        'saveResults',saveResults,'adjustmentLength',adjustmentLength,...
        'noiseScaleFactor',noiseScaleFactor,'averageSpds',averageSpds,...
        'nReversals',nReversals,'nBelowThreshold',nBelowThreshold,...
        'nObserverMatches',nObserverMatches,'thresholdScaleFactor',...
        thresholdScaleFactor);
    
    % Loop through match files and extract the number of adjustments
    dirFC = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'matchFiles',subjIDFC);
    dirThreshold = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'matchFiles',subjIDThreshold);
    countAdjustmentsFC = 0;
    countAdjustmentsThresh = 0;
    
    for j = 1:length(testWls)
        % Get directories
        simFileFC = [subjIDFC,'_',num2str(j),'.mat'];
        simFileThresh = [subjIDThreshold,'_',num2str(j),'.mat'];
        simFilePathFC = fullfile(dirFC,simFileFC);
        simFilePathThresh = fullfile(dirThreshold,simFileThresh);
        
        % Load data
        dataFC = load(simFilePathFC);
        dataThresh = load(simFilePathThresh);
        
        % Count number of adjustments
        [nTrialsFC,~] = size(dataFC.subjectSettings);
        [nTrialsThresh,~] = size(dataThresh.subjectSettings);
        
        % Add trial data to running counts 
        countAdjustmentsFC = countAdjustmentsFC + nTrialsFC;
        countAdjustmentsThresh = countAdjustmentsThresh + nTrialsThresh;
    end 
    
    nAdjustmentsFC(i) = countAdjustmentsFC; 
    nAdjustmentsThresh(i) = countAdjustmentsThresh; 
end 

% Average over observers
meanTrialsFC = mean(nAdjustmentsFC);
meanTrialsThresh = mean(nAdjustmentsThresh);