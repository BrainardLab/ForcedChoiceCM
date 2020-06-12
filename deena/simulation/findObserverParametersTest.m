% Various intial tests of using findObserverParameters to identify
% observer cone paramters.

%% Case 1 - standard observer
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_7.mat';
[testSpd, primarySpd] = getSingleMatchData(fName);
[standardFC, standardFCErr] = findObserverParameters(testSpd,primarySpd); 

% Adjustment Rule 
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_8.mat';
[testSpd, primarySpd] = getSingleMatchData(fName);
[standardAdjust, standardAdjustErr] = findObserverParameters(testSpd,primarySpd);

%% Case 2 - +2L 
% Forced choice rule 
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_9.mat';
[testSpd, primarySpd] = getSingleMatchData(fName);
[L, LErr] = findObserverParameters(testSpd,primarySpd); 

% Adjustment Rule 

%% Case 3 - +2L, -2M
% Forced choice rule 
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_10.mat';
[testSpd, primarySpd] = getSingleMatchData(fName);
[LM, LMErr] = findObserverParameters(testSpd,primarySpd); 

% Adjustment Rule 

%% Case 4 - 20% OD increase for M cone 
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_11.mat';
[testSpd, primarySpd] = getSingleMatchData(fName);
[MOD, MODErr] = findObserverParameters(testSpd,primarySpd); 

% Adjustment Rule  

%% Case 5 - 20% OD decrease for L cone, +2nm for L cone 
% Forced choice rule 
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_12.mat';
[testSpd, primarySpd] = getSingleMatchData(fName);
[LBoth, LBothErr] = findObserverParameters(testSpd,primarySpd); 

% Adjustment Rule
save('paramSearchTest.mat'); 