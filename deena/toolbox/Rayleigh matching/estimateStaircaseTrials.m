% This scripts estimates the number of trials used in a staircase Rayleigh
% matching procedure (using getMatchSeries)

% History 
%   11/2/20  dce   Wrote it.

% Rayleigh match settings (can be modified).
nObservers = 20;
coneParamsToVary = [0 0 1 1 0 1 1 0];
baseConeParams = zeros(1,length(coneParamsToVary));
opponentParams = [40.3908 205.7353 62.9590 1.0000];
age = 32;
fieldSize = 2;
rayleighPlots = false;
saveResults = false;
adjustmentLength = 51;
noiseScaleFactor = 0;  % Start with the noiseless case
monochromatic = true;
nObserverMatches = 1;
averageSpds = true;
nReversals = [1 4];
nBelowThreshold = 1;
thresholdScaleFactor = 0.5;
baseSubjID = 'TestNTrialsFiner2';
lambdaRef = 0.8;

% 8 SD matrix
  limMatrix = [570.0000   0    0.1264    0.0399    0.0459;...
  590.0000    0.0456    0.4746    0.0462    0.0716;...
  610.0000    0.2617    0.8120    0.0695    0.1325;...
  630.0000    0.6046    0.9604    0.1619    0.2685;...
  650.0000    0.8688    0.9938    0.5109    0.6458];

% 10 SD matrix 
% limMatrix =   [570.0000   0    0.1391    0.0394    0.0463;...
%   590.0000    0.0039    0.5137    0.0442    0.0733;...
%   610.0000    0.2079    0.8633    0.0647    0.1366;...
%   630.0000    0.5705    0.9930    0.1534    0.2759;...
%   650.0000    0.8574    1.0046    0.5002    0.6550];

% 6 SD matrix
% limMatrix = [570.0000,0.0090,0.1128,0.0410,0.0449;...
%              590.0000,0.0866,0.4278,0.0510,0.0673;...
%              610.0000,0.3140,0.7525,0.0809,0.1223;...
%              630.0000,0.6397,0.9215,0.1804,0.2523;...
%              650.0000,0.8814,0.9799,0.5329,0.6267];
       

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
        thresholdScaleFactor,'lambdaRef',lambdaRef,'stimLimits',limMatrix);
    
%     getMatchSeries(subjIDThreshold,sampledConeParams(i,:),opponentParams,...
%         p1Wl,p2Wl,testWls,'threshold','age',age,'fieldSize',fieldSize,...
%         'monochromatic',monochromatic,'p1Scale',p1Scale,'p2Scale',p2Scale,...
%         'testScale',testScale,'sPredicted',S,'rayleighPlots',rayleighPlots,...
%         'saveResults',saveResults,'adjustmentLength',adjustmentLength,...
%         'noiseScaleFactor',noiseScaleFactor,'averageSpds',averageSpds,...
%         'nReversals',nReversals,'nBelowThreshold',nBelowThreshold,...
%         'nObserverMatches',nObserverMatches,'thresholdScaleFactor',...
%         thresholdScaleFactor);
    
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
%         dataThresh = load(simFilePathThresh);
        
%         Count number of adjustments
        [nTrialsFC,~] = size(dataFC.subjectSettings);
%         [nTrialsThresh,~] = size(dataThresh.subjectSettings);
        
%         Add trial data to running counts 
        countAdjustmentsFC = countAdjustmentsFC + nTrialsFC;
%         countAdjustmentsThresh = countAdjustmentsThresh + nTrialsThresh;
    end 
    
    nAdjustmentsFC(i) = countAdjustmentsFC; 
%     nAdjustmentsThresh(i) = countAdjustmentsThresh; 
end 

% Average over observers
meanTrialsFC = mean(nAdjustmentsFC);
% meanTrialsThresh = mean(nAdjustmentsThresh);