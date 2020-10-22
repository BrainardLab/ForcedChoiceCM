function qpRayleighSim()
% QUEST+ Rayleigh matching experiment 
% Description:
%    This function uses the QUEST+ routine to conduct a Rayleigh matching 
%    experiment and use it to recover observer individual difference 
%    parameters. 

% History 
%    10/13/20   dce   -Wrote it 
%    10/21/20   dce   -Changed to use a single QUEST object

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

testWls = [600]; % Possible wavelengths for the test light
testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1;
end

% Set up the simulated observer. 
% coneVec is a vector of eight individual difference parameters: 
% lens pigment density, macular pigment density, L/M/S photopigment
% densities, and L/M/S peak spectral sensitivities (lambda max)

% coneVec = sampleRayleighObservers(1,zeros(1,8),[0 0 0 0 0 1 0 0]);
coneVec = [0 0 0 0 0 -3 2 0];
observer = genRayleighObserver('coneVec',coneVec);
opponentVec = [observer.colorDiffParams.lumWeight,...
    observer.colorDiffParams.rgWeight,observer.colorDiffParams.byWeight,...
    observer.colorDiffParams.noiseSd];

% Inputs
lambdaRef = 0.8;
indDiffSds = [18.7 36.5 9.0 9.0 7.4 2.0 1.5 1.3];

stimParamsDomainList = {0:0.05:1,0:0.05:1,testWls};
% psiParamsList = {0,0,-3:0.2:3,-3:0.2:3,0,-3:0.2:3,-3:0.2:3,0};
psiParamsDomainList = {0,0,0,0,0,-3:0.2:3,-3:0.2:3,0};
for i = 1:length(indDiffSds)
    psiParamsDomainList{i} = psiParamsDomainList{i}*indDiffSds(i);
end

% Psychometric function
noiseScaleFactor = 5;
PFSim = @(stimParams,coneParams)qpPFRMFull(stimParams,coneParams,...
    opponentVec,observer.colorDiffParams.noiseSd*noiseScaleFactor,S,...
    p1Spd,p2Spd,testSpds,testWls,lambdaRef);

% Simulated observer function, as a function of stimulus parameter
simObserverFun = @(stimParams) qpSimulatedObserver(stimParams,PFSim,coneVec);

% Set up a Quest object, with an option to use precomputed data if
% available
USE_PRECOMPUTE = false;
if (~USE_PRECOMPUTE)
    startTime = tic;
    fprintf('Initializing quest structure ...\n');
    questDataRaw = qpInitialize('nOutcomes', 4, ...
        'qpPF',PFSim,...
        'stimParamsDomainList',stimParamsDomainList, ...
        'psiParamsDomainList',psiParamsDomainList, ...
        'verbose', true);
    elapsedTime = toc(startTime);
    stimParamsDomainListCheck = stimParamsDomainList;
    psiParamsDomainListCheck = psiParamsDomainList;
    fprintf('Done initializing in %0.1f seconds\n',elapsedTime);
    save('questDataRaw','questDataRaw','stimParamsDomainListCheck','psiParamsDomainListCheck','-v7.3');
else
    fprintf('Loading quest structure ...\n');
    load questDataRaw questDataRaw stimParamsDomainListCheck psiParamsDomainListCheck;
    if (length(stimParamsDomainList) ~= length(stimParamsDomainListCheck))
        error('Change in stim parameters since cache of questDataRaw');
    end
    for ii = 1:length(stimParamsDomainList)
        if (any(stimParamsDomainList{ii} ~= stimParamsDomainListCheck{ii}))
            error('Change in stim parameters since cache of questDataRaw');
        end
    end
    if (length(psiParamsDomainList) ~= length(psiParamsDomainListCheck))
        error('Change in psi parameters since cache of questDataRaw');
    end
    for ii = 1:length(psiParamsDomainList)
        if (any(psiParamsDomainList{ii} ~= psiParamsDomainListCheck{ii}))
            error('Change in psi parameters since cache of questDataRaw');
        end
    end
    fprintf('Done loading\n');
end
questData = questDataRaw;

% Run trials
nTrials = 128;
for tt = 1:nTrials
    % Get stimulus for this trial
    stim = qpQuery(questData);
    
    % Simulate outcome
    outcome = simObserverFun(stim);
    
    % Update quest data structure
    questData = qpUpdate(questData,stim,outcome);
end
disp('Done Looping');
save('qDataFull','questData');

%% To do - Maximum likelihood function to estimate parameters
end