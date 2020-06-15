function [testSpds,primarySpds] = ...
    OLRayleighMatchSimSeries(subjID,observerParams,useAdjustmentRule,p1Wls,p2Wls,testWls, varargin)
% Runs a series of Rayleigh match simulations with various primary/test
% lights

% Input parsing
p = inputParser;
p.addParameter('foveal',true,@(x)(islogical(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('numReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',2,@(x) (isnumeric(x)));
p.parse(varargin{:});

% Set up subject directory 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
else 
%     error('Subject directory already exists'); 
end

% Generate an array of wavelength combinations - first column is p1, 
% second is p2, third is test
lightCombos = combvec(p1Wls,p2Wls,testWls)'; 

% Run a series of Rayleigh match experiments 
nSims = length(p1Wls)*length(p2Wls)*length(testWls);
for i = 1:nSims
    OLRayleighMatch(subjID,i,'simObserver',true,'thresholdMatching',...
        useAdjustmentRule,'observerParams',observerParams,'foveal',...
        p.Results.foveal,'p1',lightCombos(i,1),'p2',lightCombos(i,2),...
        'test',lightCombos(i,3),'nObserverMatches',...
        p.Results.nObserverMatches,'numReversals',p.Results.numReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,...
        'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
        'p2Scale',0.3,'testScale',0.5); 
end 
    
% Extract the test and match spds for each simulation
testSpds = []; 
primarySpds = []; 
for j = 1:nSims
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
    fileName = [subjID,'_',num2str(j),'.mat'];
    filePath = fullfile(outputDir,fileName);
    [testSpd,primarySpd] = getSingleMatchData(filePath); 
    testSpds = [testSpds,testSpd];
    primarySpds = [primarySpds,primarySpd]; 
end 
end 