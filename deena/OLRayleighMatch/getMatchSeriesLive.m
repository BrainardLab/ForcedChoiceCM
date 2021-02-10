function [estObserverParams, testSpds,primarySpds,testIntensities,primaryRatios] = ...
    getMatchSeriesLive(subjID,p1Wls,p2Wls,testWls,varargin)
% Runs a psychophysical experiment with a series of Rayleigh matches. 
% Syntax:
%   getMatchSeriesLive(subjID,p1Wls,p2Wls,testWls)
%
% Description:
%    Takes in information about desired primary and test wavelengths, and
%    runs OneLight Rayleigh matching experiments for the given observer
%    with each possible combination of the specified primary/test
%    wavelengths. Each set of lights is tested multiple times (as specified
%    by nObserverMatches, and pairs of lights are presented in random order. 
%    The program uses subject match spds to estimate cone individual 
%    difference parameters. 
%
% Inputs:
%    subjID             -Subject ID, entered as a character vector
%    p1Wls              -Integer or numeric vector of desired wavelengths
%                        for the first primary.
%    p2Wls              -Integer or numeric vector of desired wavelengths
%                        for the second primary.
%    testWls            -Integer or numeric vector of desired wavelengths
%                        for the test light.
%
% Outputs:
%
%
% Optional key-value pairs:
%    'fieldSize'         -Integer for observer field size. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'opponentParams'    -1x4 vector of opponent contrast parameters.  
%                         Default is [40.3908  205.7353   62.9590  1.0000].
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.2.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.2.
%    'adjustmentLength'  -Integer defining the size of the lights array
%                         available for OLRayleighMatch. Default is 201.
%    'nObserverMatches'  -Number of matches to simulate for each set of
%                         lights. Default is 1.
%    'averageSpds'       -Logical indicating to return the average of all
%                         spds for a given test wavelength. Default is
%                         true.
%    'nReversals'        -Number of reversals the observer must make before
%                         changing step size. Enter as a 2-element vector -
%                         the first element is the number of reversals for
%                         intermediate step sizes, the second is the number
%                         needed for the smallest step size. Default is
%                         [1 4].
%    'rayleighPlots'     -Logical indicating to make OLRayleighMatch plots.
%                         Default is true.
%    'saveResults'       -Logical indicating to save a file with the series
%                         of match spds. Default is true.
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on 
%                         observer parameters. Each row represents a given  
%                         test wavelength and the limits which are associated 
%                         with it. The columns are arranged as follows: 
%                         [test wl, min lambda, max lambda, min test 
%                         intensity, max test intensity]. Default is [].
%                         Note that these limits are not applied to nominal
%                         or best available matches. 
%    'LMEqualOD'         -Logical. If true, the parameter search constrains
%                         L and M cone optical densities to be equal.
%                         Default is false.
%    'dlens0'            -Logical. If true, the parameter search constrains
%                         the lens pigment density to 0. Default is true.
%    'dmac0'             -Logical. If true, the parameter search constrains
%                         macular pigment density to 0. Default is true.
%    'restrictBySd'      -Logical. If true, the parameter search restricts
%                         all params to within three standard deviations of
%                         their means. Default is true.
%    'matchErrScalar'    -Numeric scale factor to improve parameter search.
%                         Default is 100.
%    'resetAnnulus'      -logical indicating to run a script that lets the
%                         experimenter reset the annulus before the first 
%                         trial. Default is false.
% History:
%   02/10/20   dce   - Wrote it, adapted from getMatchSeries and
%                      sampleRayleighMatchSeries

% Input parsing
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.2,@(x)(isnumeric(x)));
p.addParameter('testScale',0.2,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',201,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('rayleighPlots',true,@(x)(islogical(x)));
p.addParameter('averageSpds',true,@(x)(islogical(x)));
p.addParameter('saveResults',true,@(x)(islogical(x)));
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('matchErrScalar',100,@(x)(isnumeric(x)));  
p.addParameter('resetAnnulus',false,@(x)(islogical(x)));
p.parse(varargin{:});

% Check that appropriate parameter limits have been entered
if ~isempty(p.Results.stimLimits)
    if any(p.Results.stimLimits(:,1) ~= testWls')
        error('Passed parameter limits do not match provided test wavelengths');
    end 
end 

% Set up subject directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end
if p.Results.saveResults
    fName = [subjID, '_allSpds.mat'];
    file = fullfile(outputDir,fName);
    if exist(file,'file')
        error('Specified file already exists');
    end
end

% Generate an array of wavelength combinations - first column is p1,
% second is p2, third is test
nCombos = length(p1Wls)*length(p2Wls)*length(testWls)*p.Results.nObserverMatches;
lightCombos = CombVec(p1Wls,p2Wls,testWls)';
lightCombosFull = repmat(lightCombos,p.Results.nObserverMatches,1);
randIndices = randperm(nCombos);
lightCombosFull = lightCombosFull(randIndices,:);

% Calculate Rayleigh matches for each of the light combinations
testSpds = [];
primarySpds = [];
testIntensities = [];
primaryRatios = [];

for i = 1:nCombos
    % Find parameter limits
    trialStimLimits = [];
    if ~isempty(p.Results.stimLimits)
        trialStimLimits = p.Results.stimLimits(p.Results.stimLimits(:,1)...
            ==lightCombosFull(i,3),2:5);
    end
    
    % Setup for resetting the annulus on the first trial if desired, but
    % not other trials. 
    if i==1
        resetAnnulus = p.Results.resetAnnulus;
    else
        resetAnnulus = false;
    end 
    
    % Run test 
    OLRayleighMatch(subjID,i,'simObserver',false,'foveal',...
        (p.Results.fieldSize<=2),'p1',lightCombosFull(i,1),'p2',...
        lightCombosFull(i,2),'test',lightCombosFull(i,3),'age',p.Results.age,...
        'nObserverMatches',1,'nReversals',p.Results.nReversals,'p2Scale',...
        p.Results.p2Scale,'testScale',p.Results.testScale,'p1Scale',...
        p.Results.p1Scale,'plotResponses',p.Results.rayleighPlots,...
        'adjustmentLength',p.Results.adjustmentLength,'opponentParams',...
        p.Results.opponentParams,'stimLimits',trialStimLimits',...
        'resetAnnulus',resetAnnulus);
  
        % Extract spds from the data file
        dataFile = [subjID,'_',num2str(i),'.mat'];
        dataFilePath = fullfile(outputDir,dataFile);
        [testSpd,primarySpd,testIntensity,primaryRatio] =...
            getMatchData(dataFilePath,'averageSpds',false);
         
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
% Estimate cone parameters based on data. If we want maximum likelihood
% estimation or averaging, these would be added here. 
estObserverParams = findObserverParameters(testSpds,primarySpds,...
    'age',p.Results.age,'fieldSize',p.Results.fieldSize,...
    'restrictBySd',p.Results.restrictBySd,'dlens0',p.Results.dlens0,...
    'LMEqualOD',p.Results.LMEqualOD,'dmac0',p.Results.dmac0,....
    'initialConeParams',zeros(1,8),'errScalar',p.Results.matchErrScalar,...
    'opponentParams',p.Results.opponentParams,'S',[380 2 201]);

end