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
%    'p1Scale'           -Array of numerical scale factors for the first 
%                         primary light, between 0 and 1. A different scale 
%                         factor is entered for each test wavelength. 
%                         Default is 1. Length must equal the number of
%                         test wavelengths.
%    'p2Scale'           -Array of numerical scale factors for the second 
%                         primary light, between 0 and 1. A different scale 
%                         factor is entered for each test wavelength. 
%                         Default is 0.02. Length must equal the number of
%                         test wavelengths.
%    'testScale'         -Array of numerical scale factors for the test 
%                         light, between 0 and 1. A different scale 
%                         factor is entered for each test wavelength. 
%                         Default is 0.1. Length must equal the number of
%                         test wavelengths.
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
%    'rayleighPlots'     -Logical indicating to make (and save) 
%                         OLRayleighMatch plots.Default is true.
%    'stimLimits'        -length(testWls) x 5 matrix for storing limits on
%                         observer parameters. Each row represents a given
%                         test wavelength and the limits which are associated
%                         with it. The columns are arranged as follows:
%                         [test wl, min lambda, max lambda, min test
%                         intensity, max test intensity]. Default is [].
%    'resetAnnulus'      -logical indicating to run a script that lets the
%                         experimenter reset the annulus before the first
%                         trial. Default is false.
%    'adjustment'        -Logical. When true, uses the adjustment method
%                         rather than forced choice. Default is false.
%    'pairStepSizes'        -Logical. If true, adjusts primary and test
%                            step sizes together instead of separately. 
%                            Default is false.

% History:
%   02/10/21   dce   - Wrote it, adapted from getMatchSeries and
%                      sampleRayleighMatchSeries
%   02/13/21   dce   - Edited output
%   02/25/21   dce   - Added option to have varying scale factors for
%                      different test wavelengths
%   03/09/21   dce   - Edited input to reflect changes in OLRayleighMatch
%   03/15/21   dce   - Began counterbalancing light order, added option to
%                      pair step size adjustments

% Input parsing
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.1,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',201,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('rayleighPlots',true,@(x)(islogical(x)));
p.addParameter('stimLimits',[],@(x)(isnumeric(x)));
p.addParameter('resetAnnulus',false,@(x)(islogical(x)));
p.addParameter('adjustment',false,@(x)(islogical(x)));
p.addParameter('pairStepSizes',false,@(x)(islogical(x)));
p.parse(varargin{:});

age = p.Results.age;
fieldSize = p.Results.fieldSize;
opponentParams = p.Results.opponentParams;
p1Scale = p.Results.p1Scale;
p2Scale = p.Results.p2Scale;
testScale = p.Results.testScale;
adjustmentLength = p.Results.adjustmentLength;
adjustment = p.Results.adjustment;
nObserverMatches = p.Results.nObserverMatches;

% Check that appropriate parameter limits have been entered
if ~isempty(p.Results.stimLimits)
    if any(p.Results.stimLimits(:,1) ~= testWls')
        error('Passed parameter limits do not match provided test wavelengths');
    end
end

% Check that scale factors have been entered correctly
if length(p.Results.p1Scale)~=length(testWls) ||...
        length(p.Results.p2Scale)~=length(testWls) ||...
        length(p.Results.testScale)~=length(testWls)
    error('Scale factor vectors must include the same number of elements as the number of reference wavelengths');
end 

% Generate an array of wavelength combinations - first column is p1,
% second is p2, third is test
nCombos = length(p1Wls)*length(p2Wls)*length(testWls)*p.Results.nObserverMatches;
lightCombos = CombVec(p1Wls,p2Wls,testWls)';
lightCombosFull = repmat(lightCombos,p.Results.nObserverMatches,1);
randIndices = randperm(nCombos);
lightCombosFull = lightCombosFull(randIndices,:);

% Generate scale factors for lights, where each scale factor is paired with
% a given test wavelength. Scale factors are entered into the matrix in the
% format [p1 p2 test], with one row for each row in lightCombosFull.
lightScaleFactors = zeros(nCombos,3);
for i = 1:nCombos
    lightScaleFactors(i,1) = p.Results.p1Scale(testWls==lightCombosFull(i,3));
    lightScaleFactors(i,2) = p.Results.p2Scale(testWls==lightCombosFull(i,3));
    lightScaleFactors(i,3) = p.Results.testScale(testWls==lightCombosFull(i,3));
end 

% Set up subject directory
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'matchFiles',[subjID '_' num2str(sessionNum)]);
outputFile = fullfile(outputDir,[subjID '_' num2str(sessionNum) '_summary.mat']);
fNames = cell(1,nCombos);

% Set up data arrays
testIntensities = zeros(1,nCombos);
primaryRatios = zeros(1,nCombos);

% Calculate Rayleigh matches for each of the light combinations
refFirst = round(rand(1)); % Choose randomly which light we present first on the first match 
for i = 1:nCombos
    refFirst = ~refFirst   % Light order alternates on each tria
    % Store OneLight data filename
    fNames{i} = fullfile(outputDir,[subjID,'_',num2str(sessionNum),...
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
    OLRayleighMatch([subjID '_' num2str(sessionNum)],i,'simObserver',false,'fieldSize',...
        p.Results.fieldSize,'p1',lightCombosFull(i,1),'p2',...
        lightCombosFull(i,2),'test',lightCombosFull(i,3),'age',p.Results.age,...
        'nObserverMatches',1,'nReversals',p.Results.nReversals,'p2Scale',...
        lightScaleFactors(i,2),'testScale',lightScaleFactors(i,3),'p1Scale',...
        lightScaleFactors(i,1),'plotResponses',p.Results.rayleighPlots,...
        'savePlots',p.Results.rayleighPlots,'adjustmentLength',...
        p.Results.adjustmentLength,'opponentParams',...
        p.Results.opponentParams,'stimLimits',trialStimLimits',...
        'resetAnnulus',resetAnnulus,'silent',false,...
        'adjustment',p.Results.adjustment,'testFirst',refFirst,...
        'pairStepSizes',p.Results.pairStepSizes);
    
    % Extract match position data
    [~,~,testIntensities(i),primaryRatios(i)] =...
        getMatchData(fNames{i},'nominal',false,...
        'averageSpds',false);
    
    % Save data
    save(outputFile,'fNames','testIntensities','primaryRatios','lightCombosFull',...
        'age','fieldSize','opponentParams','p1Scale','p2Scale','testScale',...
        'adjustmentLength','adjustment','nObserverMatches','testWls');
end
end