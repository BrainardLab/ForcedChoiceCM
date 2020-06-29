function [testSpds,primarySpds] = ...
    getMatchSeries(subjID,observerParams,p1Wls,p2Wls,testWls,method,varargin)
% Finds a series of simulated Rayleigh matches
% Syntax:
%   getMatchSeries(subjID,observerParams,p1Wls,p2Wls,testWls)
%
% Description:
%    Takes in observer information and cone parameters, as well as
%    information about desired primary and test wavelengths. Then, finds
%    Rayleigh matches for the given observer with each possible combination
%    of the specified primary/test wavelengths. Returns the primary and
%    test spds which were found by the observer for each match, collated
%    into two matrices. Also saves results to a file.
%
%    There are three options for calculating the Rayleigh matches in this
%    program, at varying levels of abstraction. Users can choose either
%    computePredictedRayleighMatch (analytical calculation with
%    monochromatic spds), findNominalMatch (finds the best match from a set
%    of available primary/test lights), or OLRayleighMatch (full
%    simulation). Note that for the last option, some default parameters
%    differ from those of OLRayleighMatch.
%
% Inputs:
%    subjID             -Subject ID, entered as a character vector
%    p1Wls              -Integer or numeric vector of desired wavelengths
%                        for the first primary.
%    p2Wls              -Integer or numeric vector of desired wavelengths
%                        for the second primary.
%    testWls            -Integer or numeric vector of desired wavelengths
%                        for the test light.
%    method             -Character vector for match method. Choose either
%                        'predicted', 'bestAvailable', or 'simulated'.
% Outputs:
%    testSpds           -201 x n vector representation of the predicted
%                        spds for the chosen test lights
%    primarySpds        -201 x n vector representation of the predicted
%                        spds for the chosen primary lights
%
% Optional key-value pairs:
%    'foveal'            -logical indicating whether we are making foveal
%                         matches. Default is true.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.2.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.2.
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
%    'sPredicted'        -Wavelength sampling for cone calculations when 
%                         using the predicted match method. Comes in the
%                         form [start delta nTerms]. Default is [400 1 301].

% History:
%   06/12/20  dce       Wrote it
%   06/25/20  dce       Added key-value pairs
%   06/29/20  dce       Adapted to have choice of three matching methods

% Input parsing
p = inputParser;
p.addParameter('foveal',true,@(x)(islogical(x)));
p.addParameter('saveResults',true,@(x)(islogical(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.2,@(x)(isnumeric(x)));
p.addParameter('testScale',0.2,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('thresholdMatch',true,@(x)(islogical(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('sPredicted',[400 1 301],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Check that a correct method has been entered
if ~strcmp(method,'predicted') && ~strcmp(method,'bestAvailable') && ...
        ~strcmp(method,'simulation')
    error('Unrecognized Rayleigh match method');
end

% Set up subject directory
if p.Results.saveResults
    outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),subjID);
    if (~exist(outputDir,'dir'))
        mkdir(outputDir);
    end
    fName = [subjID, '_', method, '_allSpds.mat'];
    file = fullfile(outputDir,fName);
%     if exist(file,'file')
%         error('Specified file already exists');
%     end
end

% Set field size
if p.Results.foveal
    fieldSize = 2;
else
    fieldSize = 10;
end

% Generate an array of wavelength combinations - first column is p1,
% second is p2, third is test
lightCombos = combvec(p1Wls,p2Wls,testWls)';
nCombos = length(p1Wls)*length(p2Wls)*length(testWls);

% Calculate Rayleigh matches for each of the light combinations
testSpds = [];
primarySpds = [];
for i = 1:nCombos
    if strcmp(method,'simulation')  % Use OLRayleighMatch
        % Run simulation
        OLRayleighMatch(subjID,i,'simObserver',true,'thresholdMatching',...
            p.Results.thresholdMatch,'observerParams',observerParams,'foveal',...
            p.Results.foveal,'p1',lightCombos(i,1),'p2',lightCombos(i,2),...
            'test',lightCombos(i,3),'age',p.Results.age,'nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'nBelowThreshold',p.Results.nBelowThreshold,...
            'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
            'p2Scale',p.Results.p2Scale,'testScale',p.Results.testScale,...
            'p1Scale',p.Results.p1Scale);
        % Extract spds from the data file
        simFile = [subjID,'_',num2str(i),'.mat'];
        simFilePath = fullfile(outputDir,simFile);
        [testSpd,primarySpd] = getMatchData(simFilePath);
        
    elseif strcmp(method,'bestAvailable')  % Use findNominalMatch
        % What is the name of the file we're looking for?
        lightFile = sprintf('OLRayleighMatchFineSpectralSettings_%g_%g_%g_%g_%g_%g.mat',...
            lightCombos(i,1),lightCombos(i,2),lightCombos(i,3),...
            p.Results.p1Scale,p.Results.p2Scale,p.Results.testScale);
        lightFileName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
            'precomputedStartStops',lightFile);
        % If the lights file does not exist, create it
        if ~exist(lightFileName,'file')
            OLRayleighMatchLightSettings(lightCombos(i,1),lightCombos(i,2),...
                lightCombos(i,3),'p1ScaleFactor',p.Results.p1Scale,...
                'p2ScaleFactor',p.Results.p2Scale,'testScaleFactor',...
                p.Results.testScale);
        end
        % Find the ideal match among combinations of lights in the file
        [testSpd,primarySpd] = findNominalMatch(lightFileName,...
            observerParams,'age',p.Results.age,'fieldSize',fieldSize);
        
    else  % Use computePredictedRayleighMatch
        S = [400 1 301];
        [testSpd,primarySpd] = computePredictedRayleighMatch(...
            lightCombos(i,1),lightCombos(i,2),lightCombos(i,3),...
            observerParams,'age',p.Results.age,'fieldSize',fieldSize,...
            'S',S);
    end
    % Add calculated spds from the trial to the overall array, and save
    testSpds = [testSpds,testSpd];
    primarySpds = [primarySpds,primarySpd];
    if p.Results.saveResults
        save(file,'testSpds','primarySpds','lightCombos','p');
    end
end
end