% Script to demonstrate the functionality of the qpPFRM Rayleigh matching
% psychometric function. Tests variations in chromaticity, brightness, and
% test wavelength. 

%% Setup 
close all; 

% Define monochromatic spectra to use throughout the script
S = [400 1 301];
wls = SToWls(S); 

p1Spd = zeros(length(wls),1); 
p1Wl = 670; 
p1Ind = find(wls==p1Wl);
p1Spd(p1Ind) = 1; 

p2Spd = zeros(length(wls),1);
p2Wl = 540; 
p2Ind = find(wls==p2Wl); 
p2Spd(p2Ind) = 1; 

testWls = 560:10:650; % Possible wavelengths for the test light
testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testInd = find(wls==testWls(i)); 
    testSpds(testInd,i) = 1;
end 

% Define a standard observer with default opponent contrast parameters
observer = genRayleighObserver(); 
opponentParams = [observer.colorDiffParams.lumWeight,...
    observer.colorDiffParams.rgWeight,observer.colorDiffParams.byWeight,...
    observer.colorDiffParams.noiseSd];
coneParams = ObserverParamsToVec('basic',observer);

% Set up different noise levels. These are all multiples of the noise SD
% set in the cone parameters vector, which defaults to 0.02.
noiseSd = coneParams(9);    
noise = noiseSd * [1 5 10 200];
colors = 'rgbm';       % Colors for plotting the various noise levels

%% Demo 1 - lambda
% how does p(red) change as a function of lambda?
lambdas = 0:0.02:1;
testIntensity = 1;
testWl = 600;

% 2-column vector of lambdas and stimulus parameters
stimParams = [lambdas' repmat(testIntensity,length(lambdas),1),...
    repmat(testWl,length(lambdas),1)];

% Compute probabilities for the psychometric function at different observer
% noise levels 
figure(1); 
hold on; 

for i = 1:length(noise)
    proportions = qpPFRM(stimParams,coneParams(1:8),opponentParams,...
        noise(i),S,p1Spd,p2Spd,testSpds,testWls,'judgeLum',false);
    % The first column contains probabilities of the "no" response - the
    % test light is not redder than the primary mixture, i.e. the primary
    % mixture is redder. 
    pRed = proportions(:,1);
    plot(lambdas,pRed,'o-','LineWidth',2,'Color',colors(i));
end
ylim([0 1.1]); 
xlabel('Lambda'); 
ylabel('Probability of "Primary Redder" Response');
title('qpPFRM Chromaticity Test'); 
legend(split(num2str(noise))); 

%% Demo 2 - vary test intensity 
lambda = 0.5;
testIntensities = 0:0.02:1;
testWl = 580;

% 2-column vector of lambdas and stimulus parameters
stimParams = [repmat(lambda,length(testIntensities),1),testIntensities',...
    repmat(testWl,length(lambdas),1)];

% Compute probabilities for the psychometric function at different observer
% noise levels 
figure(2); 
hold on; 

for i = 1:length(noise)
    proportions = qpPFRM(stimParams,coneParams(1:8),opponentParams,...
        noise(i),S,p1Spd,p2Spd,testSpds,testWls,'judgeLum',true);
    % The second column contains probabilities of the "yes" response - the
    % test light is brighter than the primary mixture.
    pBright = proportions(:,2);
    plot(testIntensities,pBright,'o-','LineWidth',2,'Color',colors(i));
end
ylim([0 1.1]); 
xlabel('Test Intensity'); 
ylabel('Probability of "Test Brighter" Response');
title('qpPFRM Luminance Test'); 
legend(split(num2str(noise))); 

%% Demo 3 - vary test wavelength
% Compute psychometric function for chromaticity
lambda = 0.8;
testIntensity = 0.5;
testWls = 540:2:670;
testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testInd = find(wls==testWls(i)); 
    testSpds(testInd,i) = 1;
end 

% 2-column vector of lambdas and stimulus parameters
stimParams = [repmat(lambda,length(testWls),1),...
    repmat(testIntensity,length(testWls),1),testWls'];

% Compute probabilities for the psychometric function at different observer
% noise levels 
figure(3); 
hold on; 

for i = 1:length(noise)
    proportions = qpPFRM(stimParams,coneParams(1:8),opponentParams,...
        noise(i),S,p1Spd,p2Spd,testSpds,testWls,'judgeLum',false);
    % The second column contains probabilities of the "yes" response - the  
    % test light is redder than the primary light.
    pTestRed = proportions(:,2);
    plot(testWls,pTestRed,'o-','LineWidth',2,'Color',colors(i));
end
ylim([0 1.1]); 
xlabel('Test Wavelength'); 
ylabel('Probability of "Test Redder" Response');
title('qpPFRM Test Wavelength Test'); 
legend(split(num2str(noise))); 