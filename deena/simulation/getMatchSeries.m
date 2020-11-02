function [testSpds,primarySpds,testIntensities,primaryRatios] = ...
    getMatchSeries(subjID,observerParams,opponentParams,p1Wls,p2Wls,...
    testWls,method,varargin)
% Finds a series of simulated Rayleigh matches
% Syntax:
%   getMatchSeries(subjID,observerParams,p1Wls,p2Wls,testWls)
%
% Description:
%    Takes in observer information and cone/opponent contrast parameters,
%    as well as information about desired primary and test wavelengths.
%    Then, finds Rayleigh matches for the given observer with each possible
%    combination of the specified primary/test wavelengths. Returns the
%    primary and test spds which were found by the observer for each match,
%    collated into two matrices. Also saves results to a file.
%
%    There are four options for calculating the Rayleigh matches in this
%    program, at varying levels of abstraction. Users can choose either
%    computePredictedRayleighMatch (analytical calculation with
%    monochromatic spds), searchPredictedRayleighMatch (minimize opponent contrast
%    among available spds),or OLRayleighMatch (full simulation with threshold
%    and forced-choice options). Note that when using OLRayleighMatch, some
%    defaults in this function differ from those of OLRayleighMatch.
%
% Inputs:
%    subjID             -Subject ID, entered as a character vector
%    observerParams     -Eight-element numeric vector of individual
%                        difference parameters (see ObserverVecToParams for
%                        a full description).
%    opponentParams     -Four-element numeric vector of opponent contrast
%                        parameters (see getColorDiffParams for a full
%                        description).
%    p1Wls              -Integer or numeric vector of desired wavelengths
%                        for the first primary.
%    p2Wls              -Integer or numeric vector of desired wavelengths
%                        for the second primary.
%    testWls            -Integer or numeric vector of desired wavelengths
%                        for the test light.
%    method             -Character vector for match method. Choose either
%                        'predicted', 'threshold', 'forcedChoice', or
%                        'bestAvailable'.
% Outputs:
%    testSpds           -spdLength x n vector representation of the
%                        predicted spds for the chosen test lights
%    primarySpds        -spdLength x n vector representation of the
%                        predicted spds for the chosen primary lights
%    testIntensities    -1 x n vector representation of the test intensity
%                        settings for each match, relative to the scaled
%                        spds used.
%    primaryRatios      -1 x n vector representation of the primary ratio
%                        settings for each match, relative to the scaled
%                        spds used.
%
% Optional key-value pairs:
%    'fieldSize'         -Integer for observer field size. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.2.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.2.
%    'adjustmentLength'  -Integer defining the size of the lights array
%                         available for OLRayleighMatch. Default is 3201.
%    'nObserverMatches'  -Number of matches to simulate for each set of
%                         lights. Default is 1.
%    'averageSpds'      - Logical indicating to return the average of all
%                         spds for a given test wavelength. Default is
%                         true.
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
%    'noiseScaleFactor'    -Number >=0 which determines which scalar
%                           multiple of the opponent noise SD should be
%                           used as the observer noise SD. Default is 0.
%    'nominal'           -Logical indicating to run the simulation with
%                         nominal, rather than predicted, spds. Default
%                         is false.
%    'monochromatic'     -Logical indicating to use monochromatic spds
%                         insted of those generated by the OneLight.
%                         Default is false.
%    'sPredicted'        -Wavelength sampling for cone calculations when
%                         using the predicted match method. Comes in the
%                         form [start delta nTerms]. Default is [380 2 201].
%    'rayleighPlots'     -Logical indicating to make OLRayleighMatch plots.
%                         Default is true.
%    'saveResults'       -Logical indicating to save a file with the series
%                         of match spds. Default is true.
% History:
%   06/12/20  dce       Wrote it
%   06/25/20  dce       Added key-value pairs
%   06/29/20  dce       Adapted to have choice of three matching methods
%   07/06/20  dce       Changed matching methods.
%   07/07/20  dce       Added primary ratio and test intensity as outputs
%   08/05/20  dce       Added opponent contrast info
%   08/09/20  dce       Added bestAvailable option
%   10/28/20  dce       Added OL matching with monochromatic lights

% Input parsing
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.2,@(x)(isnumeric(x)));
p.addParameter('testScale',0.2,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',3201,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('noiseScaleFactor',0,@(x)(isnumeric(x)));
p.addParameter('nominal',false,@(x)(islogical(x)));
p.addParameter('sPredicted',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('monochromatic',false,@(x)(islogical(x)));
p.addParameter('rayleighPlots',true,@(x)(islogical(x)));
p.addParameter('averageSpds',true,@(x)(islogical(x)));
p.addParameter('saveResults',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Check that a correct method has been entered
if ~strcmp(method,'predicted') && ~strcmp(method,'threshold') && ...
        ~strcmp(method,'forcedChoice') && ~strcmp(method,'bestAvailable')
    error('Unrecognized Rayleigh match method');
end

% Set up subject directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
if p.Results.saveResults
    fName = [subjID, '_', method, '_allSpds.mat'];
    file = fullfile(outputDir,fName);
    if exist(file,'file')
        error('Specified file already exists');
    end
end

% Generate an observer
if strcmp(method,'predicted') || p.Results.monochromatic % Can set S to desired value
    observer = genRayleighObserver('fieldSize',p.Results.fieldSize,'age',...
        p.Results.age,'coneVec',observerParams,'opponentParams',...
        opponentParams,'S',p.Results.sPredicted);
else                          % Leave S set to the OneLight default
    observer = genRayleighObserver('fieldSize',p.Results.fieldSize,'age',...
        p.Results.age,'coneVec',observerParams,'opponentParams',...
        opponentParams);
end

% Generate an array of wavelength combinations - first column is p1,
% second is p2, third is test
lightCombos = combvec(p1Wls,p2Wls,testWls)';
nCombos = length(p1Wls)*length(p2Wls)*length(testWls);

% Calculate Rayleigh matches for each of the light combinations
testSpds = [];
primarySpds = [];
testIntensities = [];
primaryRatios = [];
for i = 1:nCombos
    if strcmp(method,'threshold')  % Use OLRayleighMatch threshold
        % Run simulation
        OLRayleighMatch(subjID,i,'simObserver',true,'thresholdMatching',...
            true,'observerParams',observerParams,'foveal',...
            (p.Results.fieldSize<=2),'p1',lightCombos(i,1),'p2',...
            lightCombos(i,2),'test',lightCombos(i,3),'age',p.Results.age,...
            'nObserverMatches',p.Results.nObserverMatches,'nReversals',...
            p.Results.nReversals,'nBelowThreshold',...
            p.Results.nBelowThreshold,'thresholdScaleFactor',...
            p.Results.thresholdScaleFactor,'p2Scale',p.Results.p2Scale,...
            'testScale',p.Results.testScale,'p1Scale',p.Results.p1Scale,...
            'simNominalLights',p.Results.nominal,'plotResponses',...
            p.Results.rayleighPlots,'adjustmentLength',...
            p.Results.adjustmentLength,'opponentParams',opponentParams,...
            'noiseScaleFactor',p.Results.noiseScaleFactor,'monochromatic',...
            p.Results.monochromatic,'monochromaticS',p.Results.sPredicted);
        % Extract spds from the data file
        simFile = [subjID,'_',num2str(i),'.mat'];
        simFilePath = fullfile(outputDir,simFile);
        [testSpd,primarySpd,testIntensity,primaryRatio] =...
            getMatchData(simFilePath,'nominal',p.Results.nominal,...
            'averageSpds',p.Results.averageSpds);
        
    elseif strcmp(method,'forcedChoice')% Use OLRayleighMatch forced choice
        OLRayleighMatch(subjID,i,'simObserver',true,'observerParams',...
            observerParams,'foveal',(p.Results.fieldSize==2),'p1',...
            lightCombos(i,1),'p2',lightCombos(i,2),'test',lightCombos(i,3),...
            'age',p.Results.age,'nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
            'p1Scale',p.Results.p1Scale,'simNominalLights',....
            p.Results.nominal,'plotResponses',p.Results.rayleighPlots,...
            'adjustmentLength',p.Results.adjustmentLength,'opponentParams',...
            opponentParams,'noiseScaleFactor',p.Results.noiseScaleFactor,...
            'monochromatic',p.Results.monochromatic,'monochromaticS',...
            p.Results.sPredicted);
        % Extract spds from the data file
        simFile = [subjID,'_',num2str(i),'.mat'];
        simFilePath = fullfile(outputDir,simFile);
        [testSpd,primarySpd,testIntensity,primaryRatio] =...
            getMatchData(simFilePath,'nominal',p.Results.nominal,...
            'averageSpds',p.Results.averageSpds);
        
    elseif strcmp(method,'bestAvailable')  % Use bestMatchSpectra
        % Find precomputed spectra, or compute if they do not exist
        if p.Results.monochromatic  % Use monochromatic lights
            p1Scales = linspace(0,1,p.Results.adjustmentLength);
            p2Scales = 1-p1Scales;
            testScales = linspace(0,1,p.Results.adjustmentLength);
            
            wls = SToWls(p.Results.sPredicted);
            pSpds = zeros(length(wls),length(p1Scales));
            pSpds(wls==lightCombos(i,1),:) = p.Results.p1Scale*p1Scales;
            pSpds(wls==lightCombos(i,2),:) = p.Results.p2Scale*p2Scales;
            tSpds = zeros(length(wls),length(testScales));
            tSpds(wls==lightCombos(i,3),:) = p.Results.testScale*testScales;
            
        else         % Use OL spectra
            lightFile = sprintf('OLRayleighMatch%gSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
                p.Results.adjustmentLength,lightCombos(i,1),lightCombos(i,2),...
                lightCombos(i,3),p.Results.p1Scale,p.Results.p2Scale,...
                p.Results.testScale);
            lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
                'precomputedStartStops',lightFile);
            if ~exist(lightFileName,'file')
                OLRayleighMatchLightSettings(lightCombos(i,1),lightCombos(i,2),...
                    lightCombos(i,3),'p1ScaleFactor',p.Results.p1Scale,...
                    'p2ScaleFactor',p.Results.p2Scale,'testScaleFactor',...
                    p.Results.testScale,'adjustmentLength',p.Results.adjustmentLength);
            end
            lightSettings = load(lightFileName);
            testScales = lightSettings.testScales;
            p1Scales = lightSettings.testScales;
            % Identify the spds we're using
            if p.Results.nominal
                pSpds = lightSettings.primarySpdsNominal;
                tSpds = lightSettings.testSpdsNominal;
            else
                pSpds = lightSettings.primarySpdsPredicted;
                tSpds = lightSettings.testSpdsPredicted;
            end
        end
        % Find the match
        [testSpd,primarySpd,testInd,primaryInd] =...
            searchPredictedRayleighMatch(tSpds,pSpds,observer);
        % Find the primary and test ratios
        testIntensity = testScales(testInd);
        primaryRatio = p1Scales(primaryInd);
        
    else        % Use computePredictedRayleighMatch
        if p.Results.monochromatic
            p1Spd = makeMonochromaticSpd(lightCombos(i,1),...
                p.Results.p1Scale,p.Results.sPredicted);
            p2Spd = makeMonochromaticSpd(lightCombos(i,2),...
                p.Results.p2Scale,p.Results.sPredicted);
            tSpd = makeMonochromaticSpd(lightCombos(i,3),...
                p.Results.testScale,p.Results.sPredicted);
            
            [testSpd,primarySpd,testIntensity,primaryRatio] =...
                computePredictedRayleighMatch(p1Spd,p2Spd,tSpd,...
                observer,'S',p.Results.sPredicted,'addDarkSpd',false);
        else
            baseFile = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
                'precomputedSpds','OLRayleighPrimary_');
            % Names of precomputed spd files
            if p.Results.nominal
                p1File = [baseFile num2str(lightCombos(i,1)) '_'...
                    num2str(p.Results.p1Scale) '_nominal.mat'];
                p2File = [baseFile num2str(lightCombos(i,2)) '_'...
                    num2str(p.Results.p2Scale) '_nominal.mat'];
                testFile = [baseFile num2str(lightCombos(i,3)) '_'...
                    num2str(p.Results.testScale) '_nominal.mat'];
            else
                p1File = [baseFile num2str(lightCombos(i,1)) '_'...
                    num2str(p.Results.p1Scale) '_predicted_subtractDark_subtract'...
                    num2str(lightCombos(i,2)),'.mat'];
                p2File = [baseFile num2str(lightCombos(i,2)) '_'...
                    num2str(p.Results.p2Scale) '_predicted_subtractDark.mat'];
                testFile = [baseFile num2str(lightCombos(i,3)) '_'...
                    num2str(p.Results.testScale) '_predicted_subtractDark.mat'];
            end
            darkFile = [baseFile '0.mat'];
            % Load spds from the precomputed files if they exist, otherwise
            % compute
            if exist(p1File,'file')
                p1f = load(p1File);
                p1Spd = p1f.spd;
            else
                p1Spd = makeOLRayleighPrimary(lightCombos(i,1),'scaleFactor',...
                    p.Results.p1Scale,'nominal',p.Results.nominal,...
                    'subtractDark',true,'subtractAroundWl',lightCombos(i,2));
            end
            if exist(p2File,'file')
                p2f = load(p2File);
                p2Spd = p2f.spd;
            else
                p2Spd = makeOLRayleighPrimary(lightCombos(i,2),'scaleFactor',...
                    p.Results.p2Scale,'nominal',p.Results.nominal,...
                    'subtractDark',true);
            end
            if exist(testFile,'file')
                testf = load(testFile);
                tSpd = testf.spd;
            else
                tSpd = makeOLRayleighPrimary(lightCombos(i,3),'scaleFactor',...
                    p.Results.testScale,'nominal',p.Results.nominal,...
                    'subtractDark',true);
            end
            if exist(darkFile,'file')
                darkf = load(darkFile);
                darkSpd = darkf.spd;
            else
                darkSpd = makeOLRayleighPrimary(0);
            end
            
            % Run the test
            [testSpd,primarySpd,testIntensity,primaryRatio] =...
                computePredictedRayleighMatch(p1Spd,p2Spd,tSpd,...
                observer,'S',p.Results.sPredicted,'addDarkSpd',true,...
                'darkSpd',darkSpd);
        end
    end
    
    % Add calculated spds from the trial to the overall array, and save
    testSpds = [testSpds,testSpd];
    primarySpds = [primarySpds,primarySpd];
    testIntensities = [testIntensities,testIntensity];
    primaryRatios = [primaryRatios,primaryRatio];
    
    if p.Results.saveResults
        save(file,'testSpds','primarySpds','testIntensities',...
            'primaryRatios','lightCombos','p');
    end
end
end