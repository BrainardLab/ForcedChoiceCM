function [coneErrStaircase,coneErrQuest,coneErrStd] = ...
    compareQuestStaircaseRayleigh(subjID,p1Wl,p2Wl,testWls,...
    baseConeParams,coneParamsToVary,nObservers,...
    nQuestTrials,method,noiseScaleFactor,varargin)
% Runs Quest+ and staircase Rayleigh matching experiments and compares cone
% recovery
%
% Syntax:
%     compareQuestStaircaseRayleigh(subjID,p1Wl,p2Wl,testWls,...
%     baseConeParams,coneParamsToVary,nObservers,...
%     nQuestTrials,method,noiseScaleFactor)
%
% Description:
%    Runs two different Rayleigh matching procedures - staircase and QUEST+
%    - for a specified number of simulated observers, then compares the
%    accuracy of cone spectral sensitivity recovery. Both procedures are
%    run using monochromatic primary and test lights of specified
%    wavelength. The program samples a set of simulated observer
%    parameters, then uses the same set of parameters for the two versions
%    of the experiment. Many other parameters are also held constant,
%    including the possible stimulus values (lambda, test intensity, and
%    test wavelength), the observer noise level, and the base spectra. 
%
% Inputs:
%    subjID            -Character vector of subject ID
%    nObservers        -Number of simulated observers to test.
%    p1Wl              -Integer value of desired wavelength for the first
%                       primary.
%    p2Wl              -Integer value of desired wavelength for the second
%                       primary.
%    testWls           -Integer or numeric vector of desired wavelengths
%                       for the test light.
%    baseConeParams    -Eight-element numeric vector of individual
%                       difference parameters: lens pigment density,
%                       macular pigment density, L/M/S photopigment
%                       densities, and L/M/S peak spectral sensitivities
%                       (lambda max). The values are used as means for
%                       observer sampling.
%    coneParamsToVary  -Eight-element numeric vector of ones and zeros
%                       indicating which individual difference parameters
%                       should be varied. Parameters set to 1 will be 
%                       sampled around their means, while parameters set to 
%                       0 will stay at the values specified in baseParams.
%    nObservers        -Integer number of observers to simulate. 
%    nQuestTrials      -Integer number of trials to run for the QUEST+ 
%                       version of the experiment. This should be
%                       approximately equal to the number of trials in the
%                       staircase.
%    method            -Character vector indicating which method to use for
%                       staircase matching - 'predicted', 'forcedChoice',
%                      'threshold', or 'bestAvailable'.
%    noiseScaleFactor  -Number >=0 which determines which scalar multiple
%                       of the opponent noise SD should be used as the
%                       observer noise SD.
% Outputs:
%    coneErrStaircase  -1 x nObservers vector with rms differences between 
%                       sampled and recovered cone spectral sensitivities
%                       for each observer, using the staircase procedure.
%    coneErrQuest      -1 x nObservers vector with rms differences between 
%                       sampled and recovered cone spectral sensitivities
%                       for each observer, using the QUEST+ procedure.
%    coneErrStd        -1 x nObservers vector with rms differences between 
%                       sampled and standard cone spectral sensitivities
%                       for each observer.
%
% Optional key-value pairs:
%    'precomputeQuest'   -Logical. When true, the program uses precomputed
%                         QUEST data if it is available. Default is true.
%    'fieldSize'         -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'opponentParams'    -4-element vector with opponent contrast
%                         parameters. (1) is the luminance weight, (2) is
%                         the RG weight, (3) is the BY weight, and (4) is
%                         the baseline noise standard deviation. Default is
%                         [0.8078 4.1146 1.2592 0.02].
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'nStimValues'       -Number of possible stimulus values for lambda and
%                         test intensity, spaced evenly between 0 and 1. 
%                         Default is 21.
%    'lambdaRef'         -Numerical scale factor between 0 and 1 for the
%                         value of lambda used for the reference light,
%                         which is used as a baseline for computing
%                         opponent contrasts. Default is 0.8.
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'sampledObservers'  -nObservers x 8 array of previously-sampled
%                         observer parameters, useful for when observers
%                         are used across multiple experimental conditions.
%                         Default is [].
%    'nReversals'        -Number of reversals the simulated observer must
%                         make before changing step size. Enter as a 2-
%                         element vector - the first element is the number
%                         of reversals for intermediate step sizes, the
%                         second is the number needed for the smallest step
%                         size. Default is [1 4].
%    'nObserverMatches'  -Number of matches to make with each simulated set
%                         of lights. Default is 1.
%    'nBelowThreshold'   -When using a simulated observer with
%                         threshold matching, number of pairs below
%                         threshold required before recording a match.
%                         Default is 1.
%    'thresholdScaleFactor' -When using a simulated observer with
%                            threshold matching, scale factor for matching
%                            threshold. Default is 0.5.

% History
%    10/28/20   dce   -Wrote it
%    11/4/20    dce   -Changed plotting, added option to set QUEST+
%                      parameter spacing

%% Setup
% Parse input 
p = inputParser;
p.addParameter('precomputeQuest',false,@(x)(islogical(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[0.8078 4.1146 1.2592 0.0200],@(x)(isvector(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('lambdaRef',0.8,@(x)(isnumeric(x)));
p.addParameter('S',[400 1 301],@(x)(isnumeric(x)));
p.addParameter('nStimValues',51,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Define output directory 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'Quest',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
fName = [subjID '_StaircaseQuestComparison.mat'];
outputFile = fullfile(outputDir,fName);

% Get cone parameters for the simulated observers.
sampledConeParams = sampleRayleighObservers(nObservers,baseConeParams,...
    coneParamsToVary);

%% Staircase
[~,~,~,~,~,~,recoveredStaircaseParams]...
    = sampleRayleighMatch(subjID,nObservers,baseConeParams,coneParamsToVary,...
    p.Results.opponentParams,p1Wl,p2Wl,testWls,method,'age',p.Results.age,...
    'fieldSize',p.Results.fieldSize,'monochromatic',true,...
    'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
    p.Results.testScale,'S',p.Results.S,'adjustmentLength',...
    p.Results.nStimValues,'nReversals',p.Results.nReversals,...
    'nObserverMatches',p.Results.nObserverMatches,'nBelowThreshold',...
    p.Results.nBelowThreshold,'thresholdScaleFactor',...
    p.Results.thresholdScaleFactor,'noiseScaleFactor',noiseScaleFactor,...
    'LMEqualOD',false,'dlens0',true,'dmac0',true,'restrictBySd',true,...
    'makeNoObserverPlots',true,'sampledObservers',sampledConeParams);

%% QUEST+
[~,recoveredQuestParams,fittedQuestParams] = ...
    qpRayleighSim(subjID,nObservers,nQuestTrials,baseConeParams,...
    coneParamsToVary,noiseScaleFactor,p1Wl,p2Wl,testWls,'precomputeQuest',...
    p.Results.precomputeQuest,'age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'opponentParams',p.Results.opponentParams,...
    'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
    p.Results.testScale,'lambdaRef',p.Results.lambdaRef,'S',p.Results.S,...
    'plotAll',false,'plotLast',false,'sampledObservers',sampledConeParams,...
    'nStimValues',p.Results.nStimValues);

%% Analyze results
% Calculate error associated with recovered parameters
coneErrStaircase = zeros(1,nObservers);
coneErrQuest = zeros(1,nObservers);
coneErrStd = zeros(1,nObservers);
for i = 1:nObservers
    coneErrStaircase(i) = findConeSensitivityError(sampledConeParams(i,:),...
        recoveredStaircaseParams(i,:),'age',p.Results.age,'fieldSize',...
        p.Results.fieldSize,'opponentParams',p.Results.opponentParams,'S',...
        p.Results.S);
    coneErrQuest(i) = findConeSensitivityError(sampledConeParams(i,:),...
        fittedQuestParams(i,:),'age',p.Results.age,'fieldSize',...
        p.Results.fieldSize,'opponentParams',p.Results.opponentParams,'S',...
        p.Results.S);
    coneErrStd(i) = findConeSensitivityError(sampledConeParams(i,:),...
        baseConeParams,'age',p.Results.age,'fieldSize',...
        p.Results.fieldSize,'opponentParams',p.Results.opponentParams,'S',...
        p.Results.S);
end

% Make bar graph of cone errors 
coneFig = figure(); clf;
hold on;
bar([coneErrStaircase' coneErrQuest' coneErrStd']);
title('Cone Spectral Sensitivity Error');
legend('Simulated vs Staircase Params','Simulated vs QUEST+ Params',...
    'Simulated vs Standard Params');
ylabel('Error');
xlabel('Observer');

% Save results
NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_coneErrPlot.pdf']),...
    coneFig,300);
save(outputFile,'coneErrStaircase','coneErrQuest','coneErrStd','subjID','p1Wl',...
    'p2Wl','testWls','baseConeParams','coneParamsToVary','nObservers',...
    'nQuestTrials','method','noiseScaleFactor','recoveredStaircaseParams',...
    'recoveredQuestParams','fittedQuestParams');
end