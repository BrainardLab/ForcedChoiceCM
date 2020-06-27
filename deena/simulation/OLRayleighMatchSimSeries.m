function [testSpds,primarySpds] = ...
    OLRayleighMatchSimSeries(subjID,observerParams,p1Wls,p2Wls,testWls, varargin)
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
%    Note that some default parameters for this simulation differ from the
%    OLRayleighMatch default parameters.
%
% Inputs:
%    subjID             -Subject ID, entered as a character vector
%    p1Wls              -Integer or numeric vector of desired wavelengths
%                        for the first primary.
%    p2Wls              -Integer or numeric vector of desired wavelengths
%                        for the second primary.
%    testWls            -Integer or numeric vector of desired wavelengths
%                        for the test light.
% Outputs:
%    testSpds           -201 x n vector representation of the predicted
%                        spds for the chosen test lights
%    primarySpds        -201 x n vector representation of the predicted
%                        spds for the chosen primary lights
%
%
% Optional key-value pairs:
%    'foveal'            -logical indicating whether we are making foveal
%                         matches. Default is true.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.07.
%    'nObserverMatches'  -Number of matches to simulate for each set of
%                         lights. Default is 1.
%    'thresholdMatch'    -Make matches using the threshold rule, not the
%                         forced-choice rule. Default is true.
%    'nReversals'        -Number of reversals the observer must make before
%                         changing step size. Enter as a 2-element vector -
%                         the first element is the number of reversals for
%                         intermediate step sizes, the second is the number
%                         needed for the smallest step size. Default is
%                         [1 4].
%    'nBelowThreshold'   -When using a simulated observer with
%                         threshold matching, number of pairs below
%                         threshold required before recording a match.
%                         Default is 1.
%    'thresholdScaleFactor' -When using a simulated observer with
%                            threshold matching, scale factor for matching
%                            threshold. Default is 0.5.

% History:
%   06/12/20  dce       Wrote it.
%   06/25/20  dce       Added key-value pairs.

% Input parsing
p = inputParser;
p.addParameter('foveal',true,@(x)(islogical(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.07,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('thresholdMatch',true,@(x)(islogical(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x) (isnumeric(x)));
p.parse(varargin{:});

% Set up subject directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
else
    error('Subject directory already exists');
end

% Generate an array of wavelength combinations - first column is p1,
% second is p2, third is test
lightCombos = combvec(p1Wls,p2Wls,testWls)';

% Run a series of Rayleigh match experiments
nCombos = length(p1Wls)*length(p2Wls)*length(testWls);
for i = 1:nCombos
    OLRayleighMatch(subjID,i,'simObserver',true,'thresholdMatching',...
        p.Results.thresholdMatch,'observerParams',observerParams,'foveal',...
        p.Results.foveal,'p1',lightCombos(i,1),'p2',lightCombos(i,2),...
        'test',lightCombos(i,3),'nObserverMatches',...
        p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,...
        'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'p1Scale',p.Results.p1Scale);
end

% Extract the match spds from each trial
testSpds = [];
primarySpds = [];
for j = 1:nCombos
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
    fileName = [subjID,'_',num2str(i),'.mat'];
    filePath = fullfile(outputDir,fileName);
    [testSpd,primarySpd] = getMatchData(filePath);
    testSpds = [testSpds,testSpd];
    primarySpds = [primarySpds,primarySpd];
end

% Save all match spds in a single file 
file = fullfile(outputDir,'allSpds.mat');
save(file,'testSpds','primarySpds');
end