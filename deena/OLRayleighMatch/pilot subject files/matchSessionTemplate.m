%% Rayleigh matching pilot data collection session template (recreate for each subject)
% Description
%    Template script to organize runs for pilot OneLight Rayleigh match 
%    data collection. Includes code for conducting Rayleigh matches at each 
%    of six reference wavelengths in two sessions. For each reference
%    wavelength, runs two matches simultaneously - one with the primary
%    mixture displayed first, one with the reference light displayed first.
%    

% History
%    dce    xx/xx/21  - Wrote it
%    dce    06/01/21   - Changed default number of reversals
%    dce    06/02/21   - Added interleaving
%    dce    07/06/21   - Specified adjustment length and opponent params

%% To run once
% Subject information 
subjID = 'testFullProcedure2';  % Correct as needed (MELA ID)
age = 32;         % Correct as needed

% Experimental parameters
nMatchesPerSession = 1;
nReversals = [1 4];
fieldSize = 10;
p1 = 670;
p2 = 560;
isInterval = 0;
itInterval = 1;
stimInterval = 0.5;
whiteScaleFactor = 0.001;
pairStepSizes = true;
interleaveStaircases = true;
adjustmentLength = 101;
opponentParams = [40.3908  205.7353   62.9590  1.0000];

% Wavelength information 
refWls = [570 584 598 612 626 640];
p2Scalars = [0.02 0.02 0.02 0.02 0.004 0.004];
refScalars = [0.1 0.1 0.1 0.1 0.1 0.25];

% Shuffle wavelength order
rand1 = randperm(length(refWls));
rand2 = randperm(length(refWls));
shuffledWls1 = refWls(rand1);
shuffledWls2 = refWls(rand2);
shuffledP2Scalars1 = p2Scalars(rand1);
shuffledP2Scalars2 = p2Scalars(rand2);
shuffledRefScalars1 = refScalars(rand1);
shuffledRefScalars2 = refScalars(rand2);

% Save settings data 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjID);
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end 
outputFile = fullfile(outputDir,[subjID '_expSettings.mat']);
save(outputFile,'refWls','p2Scalars','refScalars','shuffledWls1','shuffledWls2',...
   'shuffledP2Scalars1','shuffledP2Scalars2','shuffledRefScalars1',...
   'shuffledRefScalars2','age','fieldSize','nMatchesPerSession','nReversals',...
   'p1','p2','isInterval','itInterval','stimInterval','whiteScaleFactor',...
   'pairStepSizes','interleaveStaircases','adjustmentLength','opponentParams');

%% To run each time you come in (before testing)
clear; 
close all;
subjID = 'testFullProcedure2';  % Correct as needed (MELA ID)

% Load settings data 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjID);
outputFile = fullfile(outputDir,[subjID '_expSettings.mat']);
settings = load(outputFile);

fNames = {};           % Data filenames
trialSessionNums = []; % Session numbers of individual trials tested in this session

%% To run each time you come in (after testing) 
% Radiometer measurements 
for i = 1:length(trialSessionNums)
    OLRadiometerMatchesPlayback(subjID,trialSessionNums(i),fNames{i},'measWhite',...
        true,'measLastOnly',false)
end 
% Close up 
ol = OneLight();
ol.shutdown;

%% Session 1 
% First ref wavelength 
trialSessionNum = 11;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls1(1),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars1(1),'testScale',settings.shuffledRefScalars1(1),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Second ref wavelength 
trialSessionNum = 12;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls1(2),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars1(2),'testScale',settings.shuffledRefScalars1(2),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Third ref wavelength 
trialSessionNum = 13;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls1(3),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars1(3),'testScale',settings.shuffledRefScalars1(3),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Fourth ref wavelength 
trialSessionNum = 14;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls1(4),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars1(4),'testScale',settings.shuffledRefScalars1(4),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Fifth ref wavelenth 
trialSessionNum = 15;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls1(5),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars1(5),'testScale',settings.shuffledRefScalars1(5),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Sixth ref wavelength 
trialSessionNum = 16;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls1(6),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars1(6),'testScale',settings.shuffledRefScalars1(6),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

%% Session 2
% First ref wavelength 
trialSessionNum = 21;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls2(1),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars2(1),'testScale',settings.shuffledRefScalars2(1),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Second ref wavelength 
trialSessionNum = 22;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls2(2),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars2(2),'testScale',settings.shuffledRefScalars2(2),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Third ref wavelength 
trialSessionNum = 23;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls2(3),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars2(3),'testScale',settings.shuffledRefScalars2(3),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Fourth ref wavelength 
trialSessionNum = 24;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls2(4),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars2(4),'testScale',settings.shuffledRefScalars2(4),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Fifth ref wavelenth 
trialSessionNum = 25;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls2(5),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars2(5),'testScale',settings.shuffledRefScalars2(5),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

% Sixth ref wavelength 
trialSessionNum = 26;  
[trialFNames, ~, ~] = getMatchSeriesLive(subjID,trialSessionNum,settings.p1,...
    settings.p2,settings.shuffledWls2(6),'nObserverMatches',...
    settings.nMatchesPerSession,'fieldSize',settings.fieldSize,...
    'adjustment',false,'nReversals',settings.nReversals,'p1Scale',1,'p2Scale',...
    settings.shuffledP2Scalars2(6),'testScale',settings.shuffledRefScalars2(6),...
    'rayleighPlots',true,'pairStepSizes',settings.pairStepSizes,'age',...
    settings.age,'isInterval',settings.isInterval,'itInterval',settings.itInterval,...
    'stimInterval',settings.stimInterval,'whiteScaleFactor',settings.whiteScaleFactor,...
    'interleaveStaircases',settings.interleaveStaircases,'adjustmentLength',...
    settings.adjustmentLength,'opponentParams',settings.opponentParams);
fNames{end+1} = trialFNames;
trialSessionNums(end+1) = trialSessionNum;

