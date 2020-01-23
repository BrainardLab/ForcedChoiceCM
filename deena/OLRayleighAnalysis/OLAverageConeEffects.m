function OLAverageConeEffects(subjectID, varargin)
% Function for calculating and plotting average cone responses across
% trials and across files

close all; 
p = inputParser;
p.addParameter('singleFile', false, @(x) (islogical(x)));
p.addParameter('measured', false, @(x) (islogical(x))); 
p.parse(varargin{:});

% Arrays for storing aggregate data. These will have three columns for L,
% M, and S cones. 
primaryCones = []; 
testCones = []; 

% Find files for subjects
baseDir= fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjectID, subjectID);
if p.Results.measured
    files = dir([baseDir, '_*_meas.mat']);
else 
   files = dir([baseDir, '_*.mat']); % how do you remove the meas files? 
end 

% Fill array with cone responses from each file 
[row, ~] = size(files); 
for i = 1:row
    fName = fullfile(files(i).folder, files(i).name);
    [pCones, tCones] = OLGetConeEffects(fName); 
    primaryCones = [primaryCones; pCones]; 
    testCones = [testCones; tCones]; 
end 

% Average primary and test cone responses 
[numMatches, ~] = size(primaryCones); 
[pLSum, pMSum, pSSum] = sum(primaryCones, 1); 
[tLSum, tMSum, tSSum] = sum(testCones, 1); 

primaryAverages = [pLSum, pMSum, pSSum] / numMatches; 
testAverages = [tLSum, tMSum, tSSum] / numMatches; 

OLPlotConeEffects(primaryAverages, testAverages, subjectID, 0, 'average', true, 'err', [tktktk]); 
end