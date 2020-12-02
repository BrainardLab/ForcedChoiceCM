% Tutorial demonstrating basic usage of the plotRayleighMatchesObserver
% program. 

% History 
%    11/18/20    dce       - Wrote it   

% Define spectra
S = [380 1 401];
wls = SToWls(S);

p1Wl = 560;
p2Wl = 670;
testWls = [570 590 610 630 650];
theLegend = {'560','590','610','630','650'};

p1Spd = zeros(length(wls),1);
p2Spd = zeros(length(wls),1);
testSpds = zeros(length(wls),length(testWls));

p1Spd(wls==p1Wl) = 1;
p2Spd(wls==p2Wl) = 1;
for i = 1:length(testWls)
    testSpds(wls==testWls(i),i) = 1;
end

% Generate standard observer
% obs = genRayleighObserver('S',S,'opponentParams',[4 2 0.5 1]);
% noise = 0.55; 
obs = genRayleighObserver('S',S);
noise = 34; 

theFig = figure();
hold on;
colors = 'rgbmy'; 
for i = 1:length(testWls)
plotRayleighMatchesObserver(obs,p1Spd,p2Spd,testSpds(:,i),noise,colors(i),...
    'Test','figHandle',theFig);
end 
legend(theLegend); 
title('Standard Observer Test');