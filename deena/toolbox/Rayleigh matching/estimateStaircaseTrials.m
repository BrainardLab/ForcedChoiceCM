function meanTrials = estimateStaircaseTrials(subjID,varargin)
% Estimates the number of trials used in a staircase Rayleigh
% matching procedure 

% Syntax:
%   estimateStaircaseTrials(subjID)
%
% Description:
%    Simulates a staircase Rayleigh matching procedure with a series of test
%    wavelengths for a specified number of observers. Computes the average
%    number of trials across observers. 
%
% Inputs:
%    subjID             -Subject ID, entered as a character vector
%
% Outputs:
%    meanTrials         -Average number of trials for the given observer
%
% Optional key-value pairs:
%    'fieldSize'         -Integer for observer field size. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'opponentParams'    -Four-element numeric vector of opponent contrast
%                         parameters (see getColorDiffParams for a full
%                         description).
%    'baseConeParams'    -Eight-element numeric vector of individual
%                         difference parameters: lens pigment density,
%                         macular pigment density, L/M/S photopigment
%                         densities, and L/M/S peak spectral sensitivities
%                         (lambda max). The values are used as means for
%                         observer sampling. Default is zeros(1,8).
%    'coneParamsToVary'  -Eight-element numeric vector of ones and zeros
%                         indicating which individual difference parameters
%                         should be varied (the noise parameter is excluded).
%                         Default is [0 0 1 1 0 1 1 0].
%    'p1Wl'              -Integer of desired wavelength for the first 
%                         primary. Default is 670.
%    'p2Wl'              -Integer of desired wavelength for the second 
%                         primary. Default is 560.
%    'testWls'           -Integer or numeric vector of desired wavelengths
%                         for the test light. Default is [570 590 610 630
%                         650].
%    'method'            -Character vector for match method. Choose either
%                         'forcedChoice' or 'threshold'.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.2.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.2.
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start delta nTerms]. Default is [380 2 201].
%    'nBelowThreshold'   -When using a simulated observer with
%                         threshold matching, number of pairs below
%                         threshold required before recording a match.
%                         Default is 1.
%    'thresholdScaleFactor' -When using a simulated observer with
%                            threshold matching, scale factor for matching
%                            threshold. Default is 0.5.
%    'noiseScaleFactor'    -Number >=0 which determines which scalar
%                           multiple of the opponent noise SD should be
%                           used as the observer noise SD. Default is 0.
%    'nObservers'        -Number of observers to sample. Default is 30.
%    'nObserverMatches'  -Number of matches to simulate for each set of
%                         lights. Default is 1.
%    'nReversals'        -Number of reversals the observer must make before
%                         changing step size. Enter as a 2-element vector -
%                         the first element is the number of reversals for
%                         intermediate step sizes, the second is the number
%                         needed for the smallest step size. Default is
%                         [1 4].
%    'monochromatic'     -Logical indicating to use monochromatic spds
%                         insted of those generated by the OneLight.
%                         Default is false.
%    'adjustmentLength'  -Integer defining the size of the lights array
%                         available for OLRayleighMatch. Default is 3201.
%    'averageSpds'       -Logical indicating to return the average of all
%                         spds for a given test wavelength. Default is
%                         true.
%    'lambdaRef'         -Number between 0 and 1 indicating which value 
%                         of lambda to use for a reference primary when
%                         calculating simulated opponent contrasts. When 
%                         empty, opponent contrasts are not computed 
%                         relative to a reference. Default is []. 
%    'limMatrix'         -length(testWls) x 5 matrix for storing limits on 
%                         stimulus parameters. Each row represents a given  
%                         test wavelength and the limits which are associated 
%                         with it. The columns are arranged as follows: 
%                         [test wl, min lambda, max lambda, min test 
%                         intensity, max test intensity]. Default is 
%                         [570.0000   0    0.1264    0.0399    0.0459;...
%                         590.0000    0.0456    0.4746    0.0462    0.0716;...
%                         610.0000    0.2617    0.8120    0.0695    0.1325;...
%                         630.0000    0.6046    0.9604    0.1619    0.2685;...
%                         650.0000    0.8688    0.9938    0.5109    0.6458].

% History 
%   11/2/20    dce   Wrote it.
%   11/21/20   dce   Turned into function  

% Parse input 
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('coneParamsToVary',[0 0 1 1 0 1 1 0],@(x)(isnumeric(x)));
p.addParameter('baseConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('p1Wl',670,@(x)(isnumeric(x)));
p.addParameter('p2Wl',560,@(x)(isnumeric(x)));
p.addParameter('testWls',[570 590 610 630 650],@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('nObservers',30,@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('noiseScaleFactor',0,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('monochromatic',true,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',101,@(x)(isnumeric(x)));
p.addParameter('averageSpds',true,@(x)(islogical(x)));
p.addParameter('lambdaRef',[],@(x)(isnumeric(x)));
p.addParameter('method','forcedChoice',@(x)(ischar(x)));
p.addParameter('limMatrix',[570.0000   0    0.1264    0.0399    0.0459;...
  590.0000    0.0456    0.4746    0.0462    0.0716;...
  610.0000    0.2617    0.8120    0.0695    0.1325;...
  630.0000    0.6046    0.9604    0.1619    0.2685;...
  650.0000    0.8688    0.9938    0.5109    0.6458],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Rayleigh match settings (can be modified).
rayleighPlots = false;
saveResults = false;

       

% Observer params
sampledConeParams = sampleRayleighObservers(p.Results.nObservers,...
    p.Results.baseConeParams,p.Results.coneParamsToVary);
nAdjustments = zeros(p.Results.nObservers,1);  % Data array

% Loop through observers and compute a series of matches
for i = 1:p.Results.nObservers
    % Compute a series of matches 
    getMatchSeries(subjID,sampledConeParams(i,:),p.Results.opponentParams,...
        p.Results.p1Wl,p.Results.p2Wl,p.Results.testWls,...
        p.Results.method,'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
        'monochromatic',p.Results.monochromatic,'p1Scale',p.Results.p1Scale,...
        'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
        'sPredicted',p.Results.S,'rayleighPlots',false,'saveResults',...
        false,'adjustmentLength',p.Results.adjustmentLength,...
        'noiseScaleFactor',p.Results.noiseScaleFactor,'averageSpds',...
        p.Results.averageSpds,'nReversals',p.Results.nReversals,...
        'nBelowThreshold',p.Results.nBelowThreshold,'nObserverMatches',...
        p.Results.nObserverMatches,'thresholdScaleFactor',...
        p.Results.thresholdScaleFactor,'lambdaRef',p.Results.lambdaRef,...
        'stimLimits',p.Results.limMatrix);
    
    % Loop through match files and extract the number of adjustments
    theDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
        'matchFiles',subjID);
    countAdjustments = 0;
    for j = 1:length(p.Results.testWls)
        % Get directories
        simFile = [subjID,'_',num2str(j),'.mat'];
        simFilePath = fullfile(theDir,simFile);
        
        % Load data
        theData = load(simFilePath);
        
        % Count number of adjustments
        [nTrials,~] = size(theData.subjectSettings);
        
        % Add trial data to running counts 
        countAdjustments = countAdjustments + nTrials;
    end 
    
    nAdjustments(i) = countAdjustments; 
end 

% Average over observers
meanTrials = mean(nAdjustments);
end 