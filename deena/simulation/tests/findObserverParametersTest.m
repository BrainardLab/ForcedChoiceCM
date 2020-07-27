% Various intial tests of using findObserverParameters to identify
% observer cone parameters.

%% Tests using OLRayleighMatch simulation and one match per observer

% Case 1 - standard observer
% Forced choice rule
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'100','100_7.mat');
[testSpd, primarySpd] = getMatchData(fName);
[standardFC, standardFCErr] = findObserverParameters(testSpd,primarySpd);

% Adjustment Rule (with default settings)
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'100','100_8.mat');
[testSpd, primarySpd] = getMatchData(fName);
[standardAdjust, standardAdjustErr] = findObserverParameters(testSpd,primarySpd);

% Case 2 - +2L
% Forced choice rule
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'100','100_9.mat');
[testSpd, primarySpd] = getMatchData(fName);
[L, LErr] = findObserverParameters(testSpd,primarySpd);

% Case 3 - +2L, -2M
% Forced choice rule
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'100','100_10.mat');
[testSpd, primarySpd] = getMatchData(fName);
[LM, LMErr] = findObserverParameters(testSpd,primarySpd);

% Case 4 - 20% OD increase for M cone
% Forced choice rule
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'100','100_11.mat');
[testSpd, primarySpd] = getMatchData(fName);
[MOD, MODErr] = findObserverParameters(testSpd,primarySpd);

% Case 5 - 20% OD decrease for L cone, +2nm for L cone
% Forced choice rule
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'100','100_12.mat');
[testSpd, primarySpd] = getMatchData(fName);
[LBoth, LBothErr] = findObserverParameters(testSpd,primarySpd);
save('paramSearchTest.mat');

%% Tests using findNominalMatch for matches
% Datasets of light settings at various test wavelengths. Search in these
% for nominal matches
% Scalars are 0.5 for test, 0.3 for p2. Test lights range from 560-640 nm
filePath = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'precomputedStartStops');
fName560 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_560_1_0.3_0.5.mat');
fName570 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_570_1_0.3_0.5.mat');
fName580 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_580_1_0.3_0.5.mat');
fName590 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_590_1_0.3_0.5.mat');
fName600 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_600_1_0.3_0.5.mat');
fName610 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_610_1_0.3_0.5.mat');
fName620 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_620_1_0.3_0.5.mat');
fName630 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_630_1_0.3_0.5.mat');
fName640 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_640_1_0.3_0.5.mat');
fName650 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_650_1_0.3_0.5.mat');
fName660 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_660_1_0.3_0.5.mat');
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
filePath = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'precomputedStartStops');
fName570 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_570_1_0.2_0.2.mat');
fName580 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_580_1_0.2_0.2.mat');
fName590 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_590_1_0.2_0.2.mat');
fName600 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_600_1_0.2_0.2.mat');
fName610 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_610_1_0.2_0.2.mat');
fName620 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_620_1_0.2_0.2.mat');
fName630 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_630_1_0.2_0.2.mat');
fName640 = fullfile(filePath,'OLRayleighMatchFineSpectralSettings_670_560_640_1_0.2_0.2.mat');
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

%% Tests using predicted matches
observerParams1 = [0 0 0 0 0 0 0 0 0];
observerParams2 = [0 0 0 0 0 4 0 0 0]; % +4nm L
observerParams3 = [0 0 0 0 0 0 -2 0 0]; % -2nm M
observerParams4 = [0 0 0 0 0 -2 2 0 0]; % shift L and M
observerParams5 = [0 0 30 0 0 0 0 0 0]; % Increase L cone OD
observerParams6 = [0 0 0 30 0 0 0 0 0]; % Increase M cone OD
observerParams7 = [0 0 20 0 0 0 0 0 0]; % Increase L cone OD
observerParams8 = [0 0 0 20 0 0 0 0 0]; % Increase M cone OD
observerParams10 = [0 0 -20 0 0 0 0 0 0]; % Change L cone OD and lambda max in opposing direction 1 
observerParams11 = [0 0 20 0 0 -2 0 0 0]; % Change L cone OD and lambda max in opposing direction 2
observerParams12 = [0 0 -20 0 0 -2 0 0 0]; % Change L cone OD and lambda max in supporting direction 1
observerParams13 = [0 0 20 0 0 2 0 0 0]; % Change L cone OD and lambda max in supporting direction 2
observerParams14 = [0 0 0 20 0 0 -2 0 0]; % Change M cone OD and lambda max in opposing direction 1
observerParams15 = [0 0 0 -20 0 0 2 0 0]; % Change M cone OD and lambda max in opposing direction 2
observerParams16 = [0 0 0 -20 0 0 -2 0 0]; % Change M cone OD and lambda max in supporting direction 1
observerParams17 = [0 0 0 20 0 0 2 0 0]; % Change M cone OD and lambda max in supporting direction 2
observerParams18 = [0 0 0 20 0 -2 0 0 0]; % Change M cone OD and L cone lambda max
observerParams = {observerParams1,observerParams2,observerParams3,...
    observerParams4,observerParams5,observerParams6,observerParams7,...
    observerParams8,observerParams10,observerParams11,...
    observerParams12,observerParams13,observerParams14,observerParams15...
    observerParams16,observerParams17,observerParams18};

calcParams = {};
S = [400 1 301];

for i = 1:length(observerParams)
    currParams = cell2mat(observerParams(i));
    initialObs = genRayleighObserver('coneVec',currParams,'S',S);
    subjID = sprintf('observerParams_%g',i);
    res = struct();
    
    [testSpds,primarySpds] = getMatchSeries(subjID,currParams,...
        670,560,570:5:660,'predicted','sPredicted',S,'saveResults',false);
    [res.observerParamsCalc,res.err] = findObserverParameters(testSpds,...
        primarySpds,'S',S);
    [res.observerParamsSetStart,res.errSetStart] = findObserverParameters(...
        testSpds,primarySpds,'initialParams',[currParams 0],'S',S);
    res.realParamsErr = findMatchError(currParams,initialObs,testSpds,...
        primarySpds,'S',S);
    res.observerParams = currParams;
    calcParams{i} = res;
end

%% Tests of predicted matches - lens/macular pigment variation 
observerParams1 = [18.7 0 0 0 0 0 0 0 0];  % dlens +1 sd
observerParams2 = [-18.7 0 0 0 0 0 0 0 0]; % dlens -1 sd
observerParams3 = [37.4 0 0 0 0 0 0 0 0];  % dlens +2 sd
observerParams4 = [-37.4 0 0 0 0 0 0 0 0]; % dlens -2 sd

observerParams5 = [0 36.5 0 0 0 0 0 0 0]; % dmac +1 sd 
observerParams6 = [0 -36.5 0 0 0 0 0 0 0]; % dmac -1 sd
observerParams7 = [0 73 0 0 0 0 0 0 0]; % dmac +2 sd
observerParams8 = [0 -73 0 0 0 0 0 0 0]; % dmac -2sd


observerParams = {observerParams1,observerParams2,observerParams3,...
    observerParams4,observerParams5,observerParams6,observerParams7,...
    observerParams8};
calcParams = {};
S = [400 1 301];

for i = 1:length(observerParams)
    currParams = cell2mat(observerParams(i));
    initialObs = genRayleighObserver('coneVec',currParams,'S',S);
    subjID = sprintf('observerParams_%g',i);
    res = struct();
    
    [testSpds,primarySpds] = getMatchSeries(subjID,currParams,...
        670,560,570:5:660,'predicted','sPredicted',S,'saveResults',false);
    [res.observerParamsCalc,res.err] = findObserverParameters(testSpds,...
        primarySpds,'S',S);
    [res.observerParamsSetStart,res.errSetStart] = findObserverParameters(...
        testSpds,primarySpds,'initialParams',[currParams 0],'S',S);
    res.realParamsErr = findMatchError(currParams,initialObs,testSpds,...
        primarySpds,'S',S);
    res.observerParams = currParams;
    calcParams{i} = res;
end

%% Tests of predicted matches - constrain cones 
observerParams1 = [0 0 20 20 0 0 0 0 0];     % Increase dphotopigment
observerParams2 = [0 0 -20 -20 0 0 0 0 0 0]; % Decrease dphotopigment
observerParams3 = [0 0 20 20 0 -2 0 0 0];    % Change dphotopigment and L cone - opposing
observerParams4 = [0 0 20 20 0 0 -2 0 0];    % Change dphotopigment and M cone - opposing
observerParams5 = [0 0 -20 -20 0 -2 0 0 0];    % Change dphotopigment and L cone - supporting
observerParams6 = [0 0 20 20 0 0 2 0 0];    % Change dphotopigment and M cone - supporting

observerParams = {observerParams1,observerParams2,observerParams3,...
    observerParams4,observerParams5,observerParams6};
calcParams = {};

for i = 1:length(observerParams)
    currParams = cell2mat(observerParams(i));
    initialObs = genRayleighObserver('coneVec',currParams,'S',[380 2 201]);
    subjID = sprintf('observerParams1_%g',i);
    res = struct();
    
    [testSpds,primarySpds] = getMatchSeries(subjID,currParams,...
        670,560,570:5:660,'predicted','sPredicted',[380 2 201]);
    [res.observerParamsCalc,res.err] = findObserverParameters(testSpds,...
        primarySpds,'LMEqualOD',true);
    [res.observerParamsSetStart,res.errSetStart] = findObserverParameters(...
        testSpds,primarySpds,'initialParams',[currParams 0],'S',S,'LMEqualOD',true);
    res.realParamsErr = findMatchError(currParams,initialObs,testSpds,...
        primarySpds);
    res.observerParams = currParams;
    calcParams{i} = res;
end

%% Tests using OLRayleighMatch for a series of matches 
observerParams1 = zeros(1,9);
observerParams2= [0 0 0 0 0 4 0 0 0]; % +4nm L
observerParams3 = [0 0 0 0 0 0 -2 0 0]; % -2nm M
observerParams4 = [0 0 0 0 0 -2 2 0 0]; % shift L and M
observerParams5 = [0 0 30 0 0 0 0 0 0]; % Increase L cone OD
observerParams6 = [0 0 -10 10 0 0 0 0 0]; % Change L and M OD

observerParams = {observerParams1,observerParams2,observerParams3,...
    observerParams4,observerParams5,observerParams6};
calcParams = {};
S = [380 2 201]; 

for i = 1:length(observerParams)
    currParams = cell2mat(observerParams(i));
    initialObs = genRayleighObserver('coneVec',currParams,'S',S);
    subjID = sprintf('observerParamsSim_%g',i);
    res = struct();
    
    [testSpds,primarySpds] = getMatchSeries(subjID,currParams,...
        670,560,570:5:650,'simulated','saveResults',false,'thresholdMatch',false);
    [res.observerParamsCalc,res.err] = findObserverParameters(testSpds,...
        primarySpds,'S',S);
    [res.observerParamsSetStart,res.errSetStart] = findObserverParameters(...
        testSpds,primarySpds,'initialParams',[currParams 0],'S');
    res.realParamsErr = findMatchError(currParams,initialObs,testSpds,...
        primarySpds,'S',S);
    res.observerParams = currParams;
    calcParams{i} = res;
end
