% Script to plot the cone responses to a subject's match settings in a 
% OneLight Rayleigh match experiment. Compares these to the cone responses 
% of the ideal match in a single bar graph. 

% History
%    dce    2/xx/20  - Wrote it
%    dce    4/5/20    - Edited for style

%% Parameters 
% Experiment data files to analyze - one with subject data, one with ideal
% matches
file = '/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_23.mat'; 
idealFile = '/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/ideal_meas.mat'; 

% Calculate subject cone responses based on measured data
measured = true; 

% Subject/session name to appear on plot 
name = 'DE_23';  

%% Compute average cone responses  
[p, t] = OLGetConeEffects(file, 'measured', measured); 
[~, numMatches] = size(p); % Number of matches in the file 

primaryAverages = mean(p, 2);
primarySD = std(p, 0, 2);

testAverages = mean(t, 2);
testSD = std(t, 0, 2);

%% Find ideal spds and their associated cone responses
% Ideal Spds (from measured data)
spd_meas = load(idealFile); 
pIdealSpd = spd_meas.primaryData(:,17); 
tIdealSpd = spd_meas.testData(:,14); 

% Cone fundamentals for standard observer  
lambdaMaxes = [558.9 530.3 420.7]';
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', 2,...
    'foveal', data.foveal); 

% Cone responses to ideal spds 
pIdealRes = T_cones * pIdealSpd; 
tIdealRes = T_cones * tIdealSpd;

%% Plot cone responses  
% Primary and test cone responses for both the subject matches and the
% ideal match 
cones = [testAverages(1), primaryAverages(1), tIdealRes(1), pIdealRes(1); ...
    testAverages(2), primaryAverages(2), tIdealRes(2), pIdealRes(2);...
    testAverages(3), primaryAverages(3), tIdealRes(3), pIdealRes(3)]; 
figure();
bar(cones);
hold on; 

% Add error bars 

% List of cones, in order they appear in the bar graph 
errCones = [cones(1,:) cones(2,:), cones(3,:)]; 
% Vector of error bars - set to 0 for ideal matches 
err = [testSD(1), primarySD(1), 0, 0, testSD(2), primarySD(2), 0, 0,...
    testSD(3), primarySD(3), 0, 0] / sqrt(numMatches);
% Axis position for error bars 
errBarPos = [0.75 0.9 1.1 1.25 1.75 1.9 2.1 2.25 2.75 2.9 3.1 3.25];
errorbar(errBarPos, errCones, err, err, 'k. ');

% Add axis labels and title 
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity');
legend('Test Average', 'Primary Average', 'Test Ideal', 'Primary Ideal');
theTitle = sprintf('%s Cone Responses vs Ideal - Measured', name);
title(theTitle); 