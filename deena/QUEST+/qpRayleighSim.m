function qpRayleighSim()
% QUEST+ Rayleigh matching experiment 
% Description:
%    This function uses the QUEST+ routine to conduct a Rayleigh matching 
%    experiment and use it to recover observer individual difference 
%    parameters. Two QUEST objects are used on alterate trials - the first 
%    tests brightness, and the second tests R/G. 

% History 
%    10/13/20   dce   -Wrote it 

% Spectra
% Define monochromatic spectra to use throughout the script
S = [400 1 301];
wls = SToWls(S);

p1Spd = zeros(length(wls),1);
p1Wl = 670;
p1Spd(wls==p1Wl) = 1;

p2Spd = zeros(length(wls),1);
p2Wl = 540;
p2Spd(wls==p2Wl) = 1;

testWls = [560 600]; % Possible wavelengths for the test light
testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1;
end

% % Set up the simulated observer
% coneVec = sampleRayleighObservers(1,zeros(1,8),[0 0 0 0 0 1 0 0]);
coneVec = [0 0 0 0 0 -3 2 0];
observer = genRayleighObserver('coneVec',coneVec);
opponentVec = [observer.colorDiffParams.lumWeight,...
    observer.colorDiffParams.rgWeight,observer.colorDiffParams.byWeight,...
    observer.colorDiffParams.noiseSd];

% Inputs
lambdaRef = 0.8;
indDiffSds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];

stimParamsDomainList = {0:0.01:1,0:0.01:1,testWls};
% psiParamsList = {0,0,-3:0.2:3,-3:0.2:3,0,-3:0.2:3,-3:0.2:3,0};
psiParamsDomainList = {0,0,-3:0.2:3,0,0,-3:0.2:3,0,0};
for i = 1:length(indDiffSds)
    psiParamsDomainList{i} = psiParamsDomainList{i}*indDiffSds(i);
end

% Psychometric functions for luminance and RG judgements
PFLum = @(stimParams,coneParams)qpPFRM(stimParams,coneParams,...
    opponentVec,observer.colorDiffParams.noiseSd*50,S,p1Spd,p2Spd,...
    testSpds,testWls,lambdaRef,'judgeLum',true);
PFRG = @(stimParams,coneParams)qpPFRM(stimParams,coneParams,...
    opponentVec,observer.colorDiffParams.noiseSd*5,S,p1Spd,p2Spd,testSpds,...
    testWls,lambdaRef,'judgeLum',false);

% Simulated observer functions, as a function of stimulus parameter
simObserverFunLum = @(stimParams) qpSimulatedObserver(stimParams,PFLum,coneVec);
simObserverFunRG = @(stimParams) qpSimulatedObserver(stimParams,PFRG,coneVec);
simObserverFuns = {simObserverFunLum simObserverFunRG};

% % Set up two Quest objects - one for luminance, one for red/green
questLumData = qpInitialize('nOutcomes', 2, ...
    'qpPF',PFLum, 'qpOutcomeF',simObserverFunLum,...
    'stimParamsDomainList',stimParamsDomainList, ...
    'psiParamsDomainList',psiParamsDomainList, ...
    'verbose', true);
questRGData = qpInitialize('nOutcomes', 2, ...
    'qpPF',PFRG,'qpOutcomeF',simObserverFunRG, ...
    'stimParamsDomainList',stimParamsDomainList, ...
    'psiParamsDomainList',psiParamsDomainList, ...
    'verbose', true);
questObjects = {questLumData questRGData};


% Run trials
nTrials = 128;
for tt = 1:nTrials
    % Index for which Quest object we're using on a given trial. Set to 1
    % on even trials, 2 on odd trials.
    qInd = mod(tt,2) + 1;
    
    % Get stimulus for this trial
    stim = qpQuery(questObjects{qInd});
    
    % Simulate outcome
    simulatedObserverFun = simObserverFuns{qInd};
    outcome = simulatedObserverFun(stim);
    
    % Update quest data structure
    questObjects{qInd} = qpUpdate(questObjects{qInd},stim,outcome);
end
disp('Done Looping');

%% To do - Maximum likelihood function to estimate parameters
end