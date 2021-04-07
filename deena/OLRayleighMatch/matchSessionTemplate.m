% Rayleigh matching pilot session template (recreate for each subject)

%% To run once
% Subject information 
subjID = 'test';  % Correct as needed
age = 32;         % Correct as needed

% Experimental parameters
nMatchesPerSession = 2; 

% Wavelength information 
refWls = [570 584 598 612 626 640];
p2Scalars = [0.02 0.02 0.02 0.02 0.0004 0.004];
refScalars = [0.1 0.1 0.1 0.1 1 1];

% Generate shuffled wavelengths
rand1 = randperm(length(refWls));
rand2 = randperm(length(refWls));
shuffledWls1 = refWls(rand1);
shuffledWls2 = refWls(rand2);

shuffledP2Scalars1 = p2Scalars(rand1);
shuffledP2Scalars2 = p2Scalars(rand2);

shuffledRefScalars1 = refScalars(rand1);
shuffledRefScalars2 = refScalars(rand2);

% Save data 

%% To run each time you come in before testing 
sessionNum = 1;   % What number time is this that you're testing?
a = [];           % Load data from before
fNames = {}; 
trialSessionNums = [];

%% To run each time you come in after testing 
% Save progress 
% Radiometer measurements 
for i = 1:length(trialSessionNums)
    OLRadiometerMatchesPlayback(subjID,trialSessionNums(i),fNames{i},'measWhite',...
        true)
end 


%% Session 1 
% First ref wavelength 
sessionNum = 01;
[trialfNames, ~, ~] = getMatchSeriesLive(subjID,sessionNum,670,560,shuffledWls1(1),...
    'nObserverMatches',nMatchesPerSession,'fieldSize',10,...
    'adjustment',false,'nReversals',[1 2],'p1Scale',1,'p2Scale',...
    shuffledP2Scalars1(1),'testScale',shuffledRefScalars1(1),'rayleighPlots',...
    true,'pairStepSizes',true,'age',age,'isInterval',0,'itInterval',1,...
    'stimInterval',0.5,'whiteScaleFactor',0.001);
% Add a new cell array of filenames
fNames{1} = trialFNames;
trialSessionNums(1) = sessionNum;

% Second ref wavelength 
% Third ref wavelength 
% Fourth ref wavelength 
% Fifth ref wavelenth 
% Sixth ref wavelength 

%% Session 2
% First ref wavelength 
% Second ref wavelength 
% Third ref wavelength 
% Fourth ref wavelength 
% Fifth ref wavelenth 
% Sixth ref wavelength
