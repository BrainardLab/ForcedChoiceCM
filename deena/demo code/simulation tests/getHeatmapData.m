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
% Start by looking at threshold data only 
% Data arrays
matchAvgErr = zeros(5,1);         % Average error across observers
matchErrs = zeros(5,20);  % Individual error of specified wls 

testWls = [610]; % Test wavelengths used in the largest increment
incs = [40 20 10 5 2];

% Data files
filePath = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'));

for i = 1:5 % 5 wl conditions 
   observersFile = load(fullfile(filePath,'paramsSearch',['varyTestWlIncrNew2_'...
       num2str(incs(i)),'_threshold'],['varyTestWlIncrNew2_'...
       num2str(incs(i)),'_threshold_paramsSearchData.mat']));
   baseObserver = observersFile.stdObs;
   errScalar = observersFile.matchErrScalar;
    for j = 1:20 % 20 observers per condition
        % Find the match spds associated with the two test wavelengths 
        primarySpds = [];
        testSpds = [];
        rayleighFile = fullfile(filePath,'matchFiles',['varyTestWlIncrNew2_'...
            num2str(incs(i)),'_threshold_' num2str(j)]);
        files = dir(rayleighFile);
        for k = 1:length(files)
            if files(k).isdir==0
               theData = load(fullfile(files(k).folder,files(k).name));
               if theData.test==testWls(1) % || theData.test==testWls(2)
                   [testSpd,primarySpd] = getMatchData(fullfile(files(k).folder,files(k).name));
                   testSpds = [testSpds,testSpd];
                   primarySpds = [primarySpds,primarySpd];
               end 
            end
        end 
        % Calculate the match error 
        observerParams = observersFile.recoveredParams(j,:);
        matchErrs(i,j) = findMatchError(observerParams,baseObserver,...
            testSpds,primarySpds,'errScalar',errScalar)/errScalar;
    end
    matchAvgErr(i) = mean(matchErrs(i,:));
end 