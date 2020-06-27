% Various intial tests of using findObserverParameters to identify
% observer cone parameters.

%% Tests using OLRayleighMatch simulation and one match per observer

% Case 1 - standard observer
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_7.mat';
[testSpd, primarySpd] = getMatchData(fName);
[standardFC, standardFCErr] = findObserverParameters(testSpd,primarySpd);

% Adjustment Rule (with default settings)
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_8.mat';
[testSpd, primarySpd] = getMatchData(fName);
[standardAdjust, standardAdjustErr] = findObserverParameters(testSpd,primarySpd);

% Case 2 - +2L
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_9.mat';
[testSpd, primarySpd] = getMatchData(fName);
[L, LErr] = findObserverParameters(testSpd,primarySpd);

% Case 3 - +2L, -2M
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_10.mat';
[testSpd, primarySpd] = getMatchData(fName);
[LM, LMErr] = findObserverParameters(testSpd,primarySpd);

% Case 4 - 20% OD increase for M cone
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_11.mat';
[testSpd, primarySpd] = getMatchData(fName);
[MOD, MODErr] = findObserverParameters(testSpd,primarySpd);

% Case 5 - 20% OD decrease for L cone, +2nm for L cone
% Forced choice rule
fName = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\100\100_12.mat';
[testSpd, primarySpd] = getMatchData(fName);
[LBoth, LBothErr] = findObserverParameters(testSpd,primarySpd);
save('paramSearchTest.mat');

%% Tests using findNominalMatch for matches
% Datasets of light settings at various test wavelengths. Search in these
% for nominal matches
% Scalars are 0.5 for test, 0.3 for p2. Test lights range from 560-640 nm
filePath = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\precomputedStartStops\';
fName560 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_560_1_0.3_0.5.mat'];
fName570 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_570_1_0.3_0.5.mat'];
fName580 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_580_1_0.3_0.5.mat'];
fName590 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_590_1_0.3_0.5.mat'];
fName600 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_600_1_0.3_0.5.mat'];
fName610 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_610_1_0.3_0.5.mat'];
fName620 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_620_1_0.3_0.5.mat'];
fName630 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_630_1_0.3_0.5.mat'];
fName640 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_640_1_0.3_0.5.mat'];
fName650 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_650_1_0.3_0.5.mat'];
fName660 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_660_1_0.3_0.5.mat'];
fNames = {fName560,fName570,fName580,fName590,fName600,fName610,fName620,...
    fName630,fName640,fName650,fName660};

% Various observer parameter vectors  (test one at a time)
observerParams = zeros(1,9); 
% observerParams = [0 0 0 0 0 4 0 0 0]; % +4nm L
% observerParams = [0 0 0 0 0 0 -2 0 0]; % -2nm M
% observerParams = [0 0 0 0 0 -2 2 0 0]; % shift L and M
% observerParams = [0 0 30 0 0 0 0 0 0]; % Increase L cone OD
% observerParams = [0 0 -10 10 0 0 0 0 0]; % Change L and M OD

% Find spds of nominal match for each of the files
testSpds = [];
primarySpds = [];
for i = 1:9
    file = (char(fNames(i)));
    [tSpd,pSpd] = findNominalMatch(file,observerParams);
    testSpds = [testSpds,tSpd];
    primarySpds = [primarySpds,pSpd];
end
checkSpdPlots = false;
if checkSpdPlots
    wls = 380:2:780;
    OLPlotSpdCheck(wls,primarySpds);
    title('Primaries');
    OLPlotSpdCheck(wls,testSpds);
    title('Test');
end
[observerParamsCalc,err] = findObserverParameters(testSpds,primarySpds);

%% Tests using findNominalMatch - take 2  
% Scalars are 0.2 for test, 0.2 for p2. Test lights range from 560-640 nm
filePath = '\Users\deena\Dropbox (Aguirre-Brainard Lab)\MELA_datadev\Experiments\ForcedChoiceCM\OLRayleighMatch\precomputedStartStops\';
fName570 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_570_1_0.2_0.2.mat'];
fName580 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_580_1_0.2_0.2.mat'];
fName590 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_590_1_0.2_0.2.mat'];
fName600 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_600_1_0.2_0.2.mat'];
fName610 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_610_1_0.2_0.2.mat'];
fName620 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_620_1_0.2_0.2.mat'];
fName630 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_630_1_0.2_0.2.mat'];
fName640 = [filePath,'OLRayleighMatchFineSpectralSettings_670_560_640_1_0.2_0.2.mat'];
fNames = {fName570,fName580,fName590,fName600,fName610,fName620,...
    fName630,fName640};
% Various observer parameter vectors  (test one at a time)
% observerParams = zeros(1,9); 
% observerParams = [0 0 0 0 0 4 0 0 0]; % +4nm L
% observerParams = [0 0 0 0 0 0 -2 0 0]; % -2nm M
% observerParams = [0 0 0 0 0 -2 2 0 0]; % shift L and M
% observerParams = [0 0 30 0 0 0 0 0 0]; % Increase L cone OD
observerParams = [0 0 -10 10 0 0 0 0 0]; % Change L and M OD

% Find spds of nominal match for each of the files
testSpds = [];
primarySpds = [];
for i = 1:length(fNames)
    file = (char(fNames(i)));
    [tSpd,pSpd] = findNominalMatch(file,observerParams);
    testSpds = [testSpds,tSpd];
    primarySpds = [primarySpds,pSpd];
end
[observerParamsCalc,err] = findObserverParameters(testSpds,primarySpds);
[observerParamsCalc2,err2] = findObserverParameters(testSpds,primarySpds,'initialParams',[observerParams 0]);