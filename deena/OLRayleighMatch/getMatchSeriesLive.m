function [fNames,testIntensities,primaryRatios] = ...
    getMatchSeriesLive(subjID,sessionNum,p1Wls,p2Wls,testWls,varargin)
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
%
% Inputs:
%    subjID             -Subject ID, entered as a character vector
%    sessionNum         -Integer session number
%    p1Wls              -Integer or numeric vector of desired wavelengths
%                        for the first primary.
%    p2Wls              -Integer or numeric vector of desired wavelengths
%                        for the second primary.
%    testWls            -Integer or numeric vector of desired wavelengths
%                        for the test light.
%
% Outputs:
%    fNames             -Cell array where each entry is the output file
%                        name for a given Rayleigh matching experiment.
%                        Used for radiometer playback. This is also saved
%                        as a file.
%    testIntensities    -1 x nCombos vector representation of the test
%                        intensity settings for each match, relative to the
%                        scaled spds used.
%    primaryRatios      -1 x nCombos vector representation of the primary
%                        ratio settings for each match, relative to the
%                        scaled spds used.
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
%    'nReversals'        -Number of reversals the observer must make before
%                         changing step size. Enter as a 2-element vector -
%                         the first element is the number of reversals for
%                         intermediate step sizes, the second is the number
%                         needed for the smallest step size. Default is
%                         [1 4].
%    'rayleighPlots'     -Logical indicating to make OLRayleighMatch plots.
%                         Default is true.
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on
%                         observer parameters. Each row represents a given
%                         test wavelength and the limits which are associated
%                         with it. The columns are arranged as follows:
%                         [test wl, min lambda, max lambda, min test
%                         intensity, max test intensity]. Default is [].
%    'resetAnnulus'      -logical indicating to run a script that lets the
%                         experimenter reset the annulus before the first
%                         trial. Default is false.
% History:
%   02/10/20   dce   - Wrote it, adapted from getMatchSeries and
%                      sampleRayleighMatchSeries
%   02/13/20   dce   - Edited output

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
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.addParameter('resetAnnulus',false,@(x)(islogical(x)));
p.parse(varargin{:});

% Check that appropriate parameter limits have been entered
if ~isempty(p.Results.stimLimits)
    if any(p.Results.stimLimits(:,1) ~= testWls')
        error('Passed parameter limits do not match provided test wavelengths');
    end
end

% Generate an array of wavelength combinations - first column is p1,
% second is p2, third is test
nCombos = length(p1Wls)*length(p2Wls)*length(testWls)*p.Results.nObserverMatches;
lightCombos = CombVec(p1Wls,p2Wls,testWls)';
lightCombosFull = repmat(lightCombos,p.Results.nObserverMatches,1);
randIndices = randperm(nCombos);
lightCombosFull = lightCombosFull(randIndices,:);

% Set up subject directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',[subjID '_' num2str(sessionNum)]);
outputFile = fullfile(outputDir,[subjID '_' num2str(sessionNum)]);
fNames = cell(1,nCombos);

% Set up data arrays
testIntensities = zeros(1,nCombos);
primaryRatios = zeros(1,nCombos);

% Calculate Rayleigh matches for each of the light combinations
for i = 1:nCombos
    % Store OneLight data filename
    fNames(i) = fullfile(outputDir,[subjID,'_',num2str(sessionNum),...
        '_',num2str(i),'.mat']);
    
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
    OLRayleighMatch([subjID '_' num2str(sessionNum)],i,'simObserver',false,'foveal',...
        (p.Results.fieldSize<=2),'p1',lightCombosFull(i,1),'p2',...
        lightCombosFull(i,2),'test',lightCombosFull(i,3),'age',p.Results.age,...
        'nObserverMatches',1,'nReversals',p.Results.nReversals,'p2Scale',...
        p.Results.p2Scale,'testScale',p.Results.testScale,'p1Scale',...
        p.Results.p1Scale,'plotResponses',p.Results.rayleighPlots,...
        'adjustmentLength',p.Results.adjustmentLength,'opponentParams',...
        p.Results.opponentParams,'stimLimits',trialStimLimits',...
        'resetAnnulus',resetAnnulus);
    
    % Extract match position data
    [~,~,testIntensities(i),primaryRatios(i)] =...
        getMatchData(fNames(i),'nominal',false,...
        'averageSpds',false);
    
    % Save data
    save(outputFile,fNames,testIntensities,primaryRatios,lightCombosFull);
    
end
end