data = load('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_19.mat');
[p, t] = OLGetConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_19.mat', 'measured', true); 
[~, numMatches] = size(p);

primaryAverages = mean(p, 2);
primarySD = std(p, 0, 2);

testAverages = mean(t, 2);
testSD = std(t, 0, 2);


% Ideal Spds
pIdealSpd = data.primarySpdsPredicted(:,198); 
tIdealSpd = data.testSpdsPredicted(:,12); 


% Cones 
lambdaMaxes = [558.9 530.3 420.7]';
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', 2,...
    'foveal', data.foveal); 

pIdealRes = T_cones * pIdealSpd; 
tIdealRes = T_cones * tIdealSpd;


figure();
cones = [primaryAverages(1), testAverages(1), pIdealRes(1), tIdealRes(1);...
    primaryAverages(2), testAverages(2), pIdealRes(2), tIdealRes(2);...
    primaryAverages(3), testAverages(3), pIdealRes(3), tIdealRes(3)]; 
bar(cones);
hold on; 

errCones = [cones(1,:) cones(2,:), cones(3,:)]; 
err = [primarySD(1), testSD(1), 0, 0, primarySD(2), testSD(2), 0, 0, primarySD(3), testSD(3), 0, 0] / sqrt(numMatches);
errBarPos = [0.75 0.9 1.1 1.25 1.75 1.9 2.1 2.25 2.75 2.9 3.1 3.25];
errorbar(errBarPos, errCones, err, err, 'k. ');


% Add axis labels
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity');
legend('Primary Average', 'Test Average', 'Primary Ideal', 'Test Ideal');
title('DE 19 Cone Responses vs Ideal - Measured');

