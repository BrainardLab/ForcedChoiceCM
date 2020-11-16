% Script for finding approximate parameter value weights for each test
% wavelength

% As of now, returns the following limit matrix for default settings:
% limMatrix = [570.0000,0.0090,0.1128,0.0410,0.0449;...
%              590.0000,0.0866,0.4278,0.0510,0.0673;...
%              610.0000,0.3140,0.7525,0.0809,0.1223;...
%              630.0000,0.6397,0.9215,0.1804,0.2523;...
%              650.0000,0.8814,0.9799,0.5329,0.6267];

% History 
%    11/12/20   dce    -Wrote it

% Observer parameters
fieldSize = 2;
age = 32;
opponentParams = [40.3908 205.7353 62.9590 1.0000];
S = [380 2 201];
nObservers = 100;
coneParamsToVary = [0 0 1 1 0 1 1 0];
baseConeParams = zeros(1,8);

% Sample 100 observers
observerParams = sampleRayleighObservers(nObservers,baseConeParams,coneParamsToVary);
observers = cell(1,nObservers);
for i = 1:nObservers
    observers{i} = genRayleighObserver('coneVec',observerParams(i,:),'age',...
        age,'fieldSize',fieldSize,'S',S,'opponentParams',opponentParams);
end 

% Wavelength parameters 
p1Wl = 670;
p2Wl = 560;
testWls = [570 590 610 630 650];
p1Scale = 1;
p2Scale = 0.02;
testScale = 0.5;

% Generate base spds
wls = SToWls(S);

p1Spd = zeros(length(wls),1);
p1Spd(wls==p1Wl) = 1*p1Scale;

p2Spd = zeros(length(wls),1);
p2Spd(wls==p2Wl) = 1*p2Scale;

testSpds = zeros(length(wls),length(testWls));
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1*testScale;
end

% Calculate nominal (analytic) matches for each observer at each test light. 
% Store nominal matches in the third dimension of the array in the format
% [lambda testIntensity];
observerMatches = zeros(length(testWls),nObservers,2);
for i = 1:length(testWls)
    for j = 1:nObservers
        [~,~,testIntensity,lambda] = ...
            computePredictedRayleighMatch(p1Spd,p2Spd,testSpds(:,i),...
            observers{j},'addDarkSpd',false,'S',S);
        observerMatches(i,j,1) = lambda;
        observerMatches(i,j,2) = testIntensity;
    end
end 

% Calculate means and standard deviations
lambdaMeans = zeros(1,length(testWls));
lambdaSDs = zeros(1,length(testWls));
tIMeans = zeros(1,length(testWls));
tISDs = zeros(1,length(testWls));

for i = 1:length(testWls)
    lambdas = observerMatches(i,:,1);
    lambdaMeans(i) = mean(lambdas);
    lambdaSDs(i) = std(lambdas);
    
    testIntensities = observerMatches(i,:,2);
    tIMeans(i) = mean(testIntensities);
    tISDs(i) = std(testIntensities);
end 

% Get ranges based on 6 standard deviations in any direction
lambdaMins = lambdaMeans-lambdaSDs*8;
lambdaMaxes = lambdaMeans+lambdaSDs*8;
tIMins = tIMeans-tISDs*8;
tIMaxes = tIMeans+tISDs*8;
limMatrix = [testWls',lambdaMins',lambdaMaxes',tIMins',tIMaxes'];

% Plot results
plotColors = 'rgbym'; % Colors for plotting
figure();
hold on;
for i = 1:length(testWls)
    for j = 1:nObservers
        % Cloud of points
        plot(observerMatches(i,j,1),observerMatches(i,j,2),'. ','Color',...
            plotColors(i));
    end 
    % 3SD rectangle
    r1 = rectangle('Position',[lambdaMeans(i)-lambdaSDs(i)*3,tIMeans(i)-tISDs(i)*3,...
        lambdaSDs(i)*6,tISDs(i)*6],'EdgeColor','k');
    
    % 6SD rectangle
    r2 = rectangle('Position',[lambdaMeans(i)-lambdaSDs(i)*6,tIMeans(i)-tISDs(i)*6,...
        lambdaSDs(i)*12,tISDs(i)*12],'EdgeColor','r');
end 
title('Distribution of Observer Matches');
xlabel('Lambda');
ylabel('Test Intensity');