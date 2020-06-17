function [testSpds,primarySpds] = ...
    OLRayleighMatchSimSeries(subjID,observerParams,useAdjustmentRule,p1Wls,p2Wls,testWls, varargin)
% Runs a series of Rayleigh match simulations with various lights
% Syntax:
%   findObserverParams(testSpds,primarySpds)
%
% Description:
%    Takes in observer information and cone parameters, as well as
%    information about desired primary and test wavelengths. Then, runs the
%    OLRayleighMatch simulation program for the given observer with each
%    possible combination of the specified primary/test wavelengths.
%    Returns the primary and test spds which were found by the observer for 
%    each match, collated into two matrices. 
%
% Inputs:
%    subjID            -Subject ID, entered as a character vector
%    useAdjustmentRule -Logical. If true, uses an adjustment matching rule.
%                       If false, uses a forced-choice matching rule.
%    p1Wls             -Integer or numeric vector of desired wavelengths
%                       for the first primary.
%    p2Wls             -Integer or numeric vector of desired wavelengths
%                       for the second primary.
%    testWls           -Integer or numeric vector of desired wavelengths
%                       for the test light.
% Outputs:
%    testSpds          -201 x n vector representation of the predicted spds 
%                       for the chosen test lights
%    primarySpds       -201 x n vector representation of the predicted spds 
%                       for the chosen primary lights
%
%
% Optional key-value pairs:
%    'foveal'          -logical indicating whether we are making foveal
%                       matches. Default is true.
%    'numReversals'   - Number of reversals the observer must make before 
%                       changing step size. Enter as a 2-element vector - 
%                       the first element is the number of reversals for 
%                       intermediate step sizes, the second is the number 
%                       needed for the smallest step size. Default is [1 4]
%    'nBelowThreshold'      - When using a simulated observer with
%                             threshold matching, number of pairs below
%                             hreshold required before recording a match.
%                             Default is 1.
%    'thresholdScaleFactor' - When using a simulated observer with
%                             threshold matching, scale factor for matching
%                             threshold. Default is 2.

% History:
%   06/12/20  dce       Wrote it.

% Input parsing
p = inputParser;
p.addParameter('foveal',true,@(x)(islogical(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('numReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',2,@(x) (isnumeric(x)));
p.parse(varargin{:});

% Set up subject directory 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
else 
%     error('Subject directory already exists'); 
end

% Generate an array of wavelength combinations - first column is p1, 
% second is p2, third is test
lightCombos = combvec(p1Wls,p2Wls,testWls)'; 

% Run a series of Rayleigh match experiments 
nSims = length(p1Wls)*length(p2Wls)*length(testWls);
for i = 1:nSims
    OLRayleighMatch(subjID,i,'simObserver',true,'thresholdMatching',...
        useAdjustmentRule,'observerParams',observerParams,'foveal',...
        p.Results.foveal,'p1',lightCombos(i,1),'p2',lightCombos(i,2),...
        'test',lightCombos(i,3),'nObserverMatches',...
        p.Results.nObserverMatches,'numReversals',p.Results.numReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,...
        'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
        'p2Scale',0.3,'testScale',0.5); 
end 
    
% Extract the test and match spds for each simulation
testSpds = []; 
primarySpds = []; 
for j = 1:nSims
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
    fileName = [subjID,'_',num2str(j),'.mat'];
    filePath = fullfile(outputDir,fileName);
    [testSpd,primarySpd] = getSingleMatchData(filePath); 
    testSpds = [testSpds,testSpd];
    primarySpds = [primarySpds,primarySpd]; 
end 
end 