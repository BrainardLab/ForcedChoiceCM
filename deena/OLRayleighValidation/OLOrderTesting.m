load('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/DE_23_ordered.mat');

% Number of trials with test and primary first
numPFirst = 5;
numTFirst = 5;

% Find means and standard errors
pAverages_pFirst = mean(p_primaryFirst, 2);
pSEM_pFirst = std(p_primaryFirst, 0, 2) / sqrt(numPFirst);

pAverages_tFirst = mean(p_testFirst, 2);
pSEM_tFirst = std(p_testFirst, 0, 2) / sqrt (numTFirst);

tAverages_pFirst = mean(t_primaryFirst, 2);
tSEM_pFirst = std(t_primaryFirst, 0, 2) / sqrt(numPFirst);

tAverages_tFirst = mean(t_testFirst, 2);
tSEM_tFirst = std(t_testFirst, 0, 2) / sqrt(numTFirst);

figure();
cones = [pAverages_pFirst(1), pAverages_tFirst(1), tAverages_pFirst(1), tAverages_tFirst(1);...
    pAverages_pFirst(2), pAverages_tFirst(2), tAverages_pFirst(2), tAverages_tFirst(2);...
    pAverages_pFirst(3), pAverages_tFirst(3), tAverages_pFirst(3), tAverages_tFirst(3)];
bar(cones);
hold on;

% Add error
errCones = [pAverages_pFirst(1), pAverages_tFirst(1), tAverages_pFirst(1), tAverages_tFirst(1),...
    pAverages_pFirst(2), pAverages_tFirst(2), tAverages_pFirst(2), tAverages_tFirst(2),...
    pAverages_pFirst(3), pAverages_tFirst(3), tAverages_pFirst(3), tAverages_tFirst(3)];

err = [pSEM_pFirst(1), pSEM_tFirst(1), tSEM_pFirst(1), tSEM_tFirst(1),...
    pSEM_pFirst(2), pSEM_tFirst(2), tSEM_pFirst(2), tSEM_tFirst(2),...
    pSEM_pFirst(3), pSEM_tFirst(3), tSEM_pFirst(3), tSEM_tFirst(3)];

errBarPos = [0.75 0.9 1.1 1.25 1.75 1.9 2.1 2.25 2.75 2.9 3.1 3.25];
errorbar(errBarPos, errCones, err, err, 'k. ');


% Add axis labels
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity');
legend('primary: primary first', 'primary: test first', 'test: primary first', 'test: test first');
title('DE 23 Measured Cone Responses: Order Effect');