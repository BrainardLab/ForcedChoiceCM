file = load('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_19.mat')
[p, t] = OLGetConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE/DE_19.mat'); 
[~, numMatches] = size(p);

primaryAverages = mean(p, 2);
% primarySD = std(p, 2);

testAverages = mean(t, 2);
% testSD = std(t, 2);

% err = [testSD(1), primarySD(1), testSD(2), primarySD(2), testSD(3), primarySD(3)] / sqrt(numMatches);


% Ideal Spds
pIdealSpd = file.primarySpdsPredicted(:,198); 
tIdealSpd = file.testSpdsPredicted(:,12); 



% Cones 
lambdaMaxes = [558.9 530.3 420.7]';
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', 2,...
    'foveal', file.foveal); 

pIdealRes = T_cones * pIdealSpd; 
tIdealRes = T_cones * tIdealSpd;


figure();
cones = [primaryAverages(1), testAverages(1), pIdealRes(1), tIdealRes(1);...
    primaryAverages(2), testAverages(2), pIdealRes(2), tIdealRes(2);...
    primaryAverages(3), testAverages(3), pIdealRes(3), tIdealRes(3)]; 
bar(cones);


% Add axis labels
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity');
legend('Primary Average', 'Test Average', 'Primary Ideal', 'Test Ideal');
title('DE 19 Cone Responses vs Ideal');

