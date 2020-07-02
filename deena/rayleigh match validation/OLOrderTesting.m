% Script to test for an order effect in the OneLight Rayleigh match 
% experiments. Takes in a session data file and groups trials based on 
% whether the primary or test light wa shown first. Then calculates and 
% plots the average primary and test cone responses for the two groups of
% trials. 

% History
%    dce    2/xx/20  - Wrote it
%    dce    4/5/20    - Edited for style

% Data file of choice 

fName = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDataDir'),DE_peripheral,DE_peripheral_1.mat');
name = 'DE_peripheral'; 

% Cone responses to chosen matches (measured) 
[pCones,tCones] = OLGetConeEffects(fName, 'measured', true);

% Find number of trials with test and primary first. Assume that the user
% began with the primary first and alternated after each match. 
[row, col] = size(pCones); 
numPFirst = ceil(col/2); 
numTFirst = floor(col/2);

% Initialize data arrays 
p_primaryFirst = zeros(3,numPFirst);
p_testFirst = zeros(3,numTFirst);
t_primaryFirst = zeros(3,numPFirst);
t_testFirst = zeros(3,numTFirst);

% Indices to track data entry 
pFirstIndex = 1;
tFirstIndex = 1;

% Loop through and calculate cone effects, place in appropriate array. Odd
% numbers are primary first, even numbers are test first 
for i = 1:col
    if mod(i,2) ~= 0 
        % Odd trials are primary first  
        p_primaryFirst(:,pFirstIndex) = pCones(:,i);
        t_primaryFirst(:,pFirstIndex) = tCones(:,i);
        pFirstIndex = pFirstIndex+1;
    else            
        % Even trials are test first 
        p_testFirst(:,tFirstIndex) = pCones(:,i);
        t_testFirst(:,tFirstIndex) = tCones(:,i);
        tFirstIndex = tFirstIndex+1; 
    end 
end 

% Find means and standard errors for the different groups 
pAverages_pFirst = mean(p_primaryFirst, 2);
pSEM_pFirst = std(p_primaryFirst, 0, 2) / sqrt(numPFirst);

pAverages_tFirst = mean(p_testFirst, 2);
pSEM_tFirst = std(p_testFirst, 0, 2) / sqrt (numTFirst);

tAverages_pFirst = mean(t_primaryFirst, 2);
tSEM_pFirst = std(t_primaryFirst, 0, 2) / sqrt(numPFirst);

tAverages_tFirst = mean(t_testFirst, 2);
tSEM_tFirst = std(t_testFirst, 0, 2) / sqrt(numTFirst);

% Make a bar plot of the averages, 
figure();
cones = [pAverages_pFirst(1), pAverages_tFirst(1), tAverages_pFirst(1), tAverages_tFirst(1);...
    pAverages_pFirst(2), pAverages_tFirst(2), tAverages_pFirst(2), tAverages_tFirst(2);...
    pAverages_pFirst(3), pAverages_tFirst(3), tAverages_pFirst(3), tAverages_tFirst(3)];
bar(cones);
hold on;

% Add error bars 
% Vector of cones in order they appear on the plot
errCones = [pAverages_pFirst(1), pAverages_tFirst(1), tAverages_pFirst(1), tAverages_tFirst(1),...
    pAverages_pFirst(2), pAverages_tFirst(2), tAverages_pFirst(2), tAverages_tFirst(2),...
    pAverages_pFirst(3), pAverages_tFirst(3), tAverages_pFirst(3), tAverages_tFirst(3)];

% Vector of error bar values in order  
err = [pSEM_pFirst(1), pSEM_tFirst(1), tSEM_pFirst(1), tSEM_tFirst(1),...
    pSEM_pFirst(2), pSEM_tFirst(2), tSEM_pFirst(2), tSEM_tFirst(2),...
    pSEM_pFirst(3), pSEM_tFirst(3), tSEM_pFirst(3), tSEM_tFirst(3)];

% Positions for error bars
errBarPos = [0.75 0.9 1.1 1.25 1.75 1.9 2.1 2.25 2.75 2.9 3.1 3.25];
errorbar(errBarPos, errCones, err, err, 'k. ');

% Add axis labels
names ={'L'; 'M'; 'S' };
set(gca,'xticklabel', names)
ylabel('Relative Response Intensity');
legend('primary: primary first', 'primary: test first', 'test: primary first', 'test: test first');
theTitle = sprintf('%s Measured Cone Responses: Order Effect', name);
title(theTitle); 