% Data arrays
coneErrFC = zeros(4);
matchErrFC = zeros(4);
coneErrAdjust = zeros(4);
matchErrAdjust = zeros(4);

% Data files
filePath = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'paramsSearch');

% Forced choice
fNameFC1 = fullfile(filePath,'testNoiseIncr2_40_0_FC\paramsSearchData.mat');
fNameFC2 = fullfile(filePath,'testNoiseIncr2_40_0.01_FC\paramsSearchData.mat');
fNameFC3 = fullfile(filePath,'testNoiseIncr2_40_0.02_FC\paramsSearchData.mat');
fNameFC4 = fullfile(filePath,'testNoiseIncr2_40_0.04_FC\paramsSearchData.mat');
fNameFC5 = fullfile(filePath,'testNoiseIncr2_20_0_FC\paramsSearchData.mat');
fNameFC6 = fullfile(filePath,'testNoiseIncr2_20_0.01_FC\paramsSearchData.mat');
fNameFC7 = fullfile(filePath,'testNoiseIncr2_20_0.02_FC\paramsSearchData.mat');
fNameFC8 = fullfile(filePath,'testNoiseIncr2_20_0.04_FC\paramsSearchData.mat');
fNameFC9 = fullfile(filePath,'testNoiseIncr2_10_0_FC\paramsSearchData.mat');
fNameFC10 = fullfile(filePath,'testNoiseIncr2_10_0.01_FC\paramsSearchData.mat');
fNameFC11 = fullfile(filePath,'testNoiseIncr2_10_0.02_FC\paramsSearchData.mat');
fNameFC12 = fullfile(filePath,'testNoiseIncr2_10_0.04_FC\paramsSearchData.mat');
fNameFC13 = fullfile(filePath,'testNoiseIncr2_5_0_FC\paramsSearchData.mat');
fNameFC14 = fullfile(filePath,'testNoiseIncr2_5_0.01_FC\paramsSearchData.mat');
fNameFC15 = fullfile(filePath,'testNoiseIncr2_5_0.02_FC\paramsSearchData.mat');
fNameFC16 = fullfile(filePath,'testNoiseIncr2_5_0.04_FC\paramsSearchData.mat');

filesFC = {fNameFC1 fNameFC2 fNameFC3 fNameFC4; fNameFC5 fNameFC6 fNameFC7...
    fNameFC8; fNameFC9 fNameFC10 fNameFC11 fNameFC12; fNameFC13 fNameFC14...
    fNameFC15 fNameFC16};
% Adjustment
fNameAdjust1 = fullfile(filePath,'testNoiseIncr2_40_0_threshold\paramsSearchData.mat');
fNameAdjust2 = fullfile(filePath,'testNoiseIncr2_40_0.01_threshold\paramsSearchData.mat');
fNameAdjust3 = fullfile(filePath,'testNoiseIncr2_40_0.02_threshold\paramsSearchData.mat');
fNameAdjust4 = fullfile(filePath,'testNoiseIncr2_40_0.04_threshold\paramsSearchData.mat');
fNameAdjust5 = fullfile(filePath,'testNoiseIncr2_20_0_threshold\paramsSearchData.mat');
fNameAdjust6 = fullfile(filePath,'testNoiseIncr2_20_0.01_threshold\paramsSearchData.mat');
fNameAdjust7 = fullfile(filePath,'testNoiseIncr2_20_0.02_threshold\paramsSearchData.mat');
fNameAdjust8 = fullfile(filePath,'testNoiseIncr2_20_0.04_threshold\paramsSearchData.mat');
fNameAdjust9 = fullfile(filePath,'testNoiseIncr2_10_0_threshold\paramsSearchData.mat');
fNameAdjust10 = fullfile(filePath,'testNoiseIncr2_10_0.01_threshold\paramsSearchData.mat');
fNameAdjust11 = fullfile(filePath,'testNoiseIncr2_10_0.02_threshold\paramsSearchData.mat');
fNameAdjust12 = fullfile(filePath,'testNoiseIncr2_10_0.04_threshold\paramsSearchData.mat');
fNameAdjust13 = fullfile(filePath,'testNoiseIncr2_5_0_threshold\paramsSearchData.mat');
fNameAdjust14 = fullfile(filePath,'testNoiseIncr2_5_0.01_threshold\paramsSearchData.mat');
fNameAdjust15 = fullfile(filePath,'testNoiseIncr2_5_0.02_threshold\paramsSearchData.mat');
fNameAdjust16 = fullfile(filePath,'testNoiseIncr2_5_0.04_threshold\paramsSearchData.mat');

filesAdjust = {fNameAdjust1 fNameAdjust2 fNameAdjust3 fNameAdjust4;...
    fNameAdjust5 fNameAdjust6 fNameAdjust7 fNameAdjust8; fNameAdjust9...
    fNameAdjust10 fNameAdjust11 fNameAdjust12; fNameAdjust13 fNameAdjust14...
    fNameAdjust15 fNameAdjust16};

for i = 1:4
    for j = 1:4
        fcFile = load(char(filesFC(i,j)));
        thresholdFile = load(char(filesAdjust(i,j)));
        
        coneErrFC(i,j) = fcFile.coneAvgErr;
        matchErrFC(i,j) = fcFile.matchAvgErr;
        coneErrAdjust(i,j) = thresholdFile.coneAvgErr;
        matchErrAdjust(i,j) = thresholdFile.matchAvgErr;
    end
end
%% test match error for wl increments 
% Data arrays
matchAvgErr = zeros(1,6);         % Average error 
matchErrs = {[],[],[],[],[],[]};  % Individual error of specified wls 

% Data files
filePath = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),'paramsSearch');
fName1 = fullfile(filePath,'varyTestWl5_1_threshold','paramsSearchData.mat');
fName2 = fullfile(filePath,'varyTestWl5_2_threshold','paramsSearchData.mat');
fName3 = fullfile(filePath,'varyTestWl5_3_threshold','paramsSearchData.mat');
fName4 = fullfile(filePath,'varyTestWl5_4_threshold','paramsSearchData.mat');
fName5 = fullfile(filePath,'varyTestWl5_5_threshold','paramsSearchData.mat');
fName6 = fullfile(filePath,'varyTestWl5_6_threshold','paramsSearchData.mat');

testWls = [570 610]; % Test wavelengths used in the largest increment

fileNames = {fName1 fName2 fName3 fName4 fName5 fName6};
for i = 1:length(fileNames)
    load(cell2mat(fileNames(i)));
end 