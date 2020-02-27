data = load('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_23.mat');
[p, t] = OLGetConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_23.mat', 'measured', true); 
[~, numMatches] = size(p);

primaryAverages = mean(p, 2);
primarySD = std(p, 0, 2);

testAverages = mean(t, 2);
testSD = std(t, 0, 2);

% Ideal Spds (from measured data
spd_meas = load('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/ideal_meas.mat')
pIdealSpd = spd_meas.primaryData(:,17); 
tIdealSpd = spd_meas.testData(:,14); 

% Cones 
lambdaMaxes = [558.9 530.3 420.7]';
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', 2,...
    'foveal', data.foveal); 

pIdealRes = T_cones * pIdealSpd; 
tIdealRes = T_cones * tIdealSpd;

figure();
cones = [testAverages(1), primaryAverages(1), tIdealRes(1), pIdealRes(1); ...
    testAverages(2), primaryAverages(2), tIdealRes(2), pIdealRes(2);...
    testAverages(3), primaryAverages(3), tIdealRes(3), pIdealRes(3)]; 
bar(cones);
hold on; 

errCones = [cones(1,:) cones(2,:), cones(3,:)]; 
err = [testSD(1), primarySD(1), 0, 0, testSD(2), primarySD(2), 0, 0,...
    testSD(3), primarySD(3), 0, 0] / sqrt(numMatches);
errBarPos = [0.75 0.9 1.1 1.25 1.75 1.9 2.1 2.25 2.75 2.9 3.1 3.25];
errorbar(errBarPos, errCones, err, err, 'k. ');

% Add axis labels
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity');
legend('Test Average', 'Primary Average', 'Test Ideal', 'Primary Ideal');
title('DE 23 Cone Responses vs Ideal - Measured');
