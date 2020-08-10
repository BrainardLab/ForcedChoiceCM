function sampleRayleighMatchSeries(subjID,nObservers,p1,p2,test,...
    coneParamsToVary,testingParamToVary,testingValsToVary,...
    varargin)
% Runs Rayleigh matching simulations with different groups of observers
% Syntax:
%   sampleRayleighMatchSeries(subjID,nObservers,p1,p2,test,coneParamsToVary,
%   testingParamToVary,testingValsToVary,varargin)
%
% Description:
%    This program tests how varying different experimental settings
%    affects the ability to recover observer parameters using
%    sampleRayleighMatch. The experimenter specifies a testing parameter to
%    systematically vary. Currently, options include observer noise 
%    (expressed as a scalar multiple of opponent noise SD), increment of 
%    test wavelengths, and number of matches per trial). For each value of 
%    the parameter, the program runs sampleRayleighMatch, which samples 
%    observers, simulates Rayleigh matching, and attempts to recover cone 
%    parameters. Both the forced-choice and threshold versions are tested
%    for each variable parameter value.
%
%    The program makes and saves plots of two types of error - average cone
%    spectral sensitivity error and average match error - allowing
%    comparison across different values of the variable parameter.
%
%    There are also options to vary two parameters at once (current options
%    include noise and test wavelength increment, number of matches and
%    test wavelength increment, and number of matches and noise. In this
%    case, comparative "heatmaps" of the two parameters are produced
%    instead of the single parameter plots.
%
% Inputs:
%    subjID         -Character vector of subject ID
%    nObservers     -Number of simulated observers to sample.
%    p1             -Integer or numeric vector of desired wavelengths
%                    for the first primary.
%    p2             -Integer or numeric vector of desired wavelengths
%                    for the second primary.
%    test           -Integer or numeric vector of desired wavelengths
%                    for the test light.
%    coneParamsToVary   -Eight-element numeric vector of ones and zeros
%                        indicating which individual difference parameters
%                        should be varied (the noise parameter is excluded).
%                        Parameters set to 1 will be sampled around their
%                        standard deviation, while parameters set to 0 will
%                        stay at the values specified in baseParams.
%    testingParamToVary -Character vector of the parameter to vary across
%                        groups of observers. Current options include
%                        'noise','nMatches', and 'testWlIncr'.
%    testingValsToVary  -Vector which includes the different values of the
%                        variable parameter to test.
%
% Outputs:
%    None (saves figures)
%
% Optional key-value pairs:
%    'fieldSize          -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'baseConeParams'    -Eight-element numeric vector of individual
%                         difference parameters used as a starting point.
%                         Default is zeros(1,8)
%    opponentParams      -4-element vector with opponent contrast  
%                         parameters. Default is [0.8078 4.1146 1.2592 0.02].
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'adjustmentLength'  -Integer defining the size of the lights array
%                         available for OLRayleighMatch. Default is 3201.
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
%    'noiseScaleFactor'  -Number >=0 which determines which scalar 
%                         multiple of the opponent noise SD should be 
%                         used as the observer noise SD. Default is 0.
%    'LMEqualOD'         -Logical. If true, the parameter search constrnMatains
%                         L and M cone optical densities to be equal.
%                         Default is false.
%    'dlens0'            -Logical. If true, the parameter search constrains
%                         the lens pigment density to 0. Default is true.
%    'dmac0'             -Logical. If true, the parameter search constrains
%                         macular pigment density to 0. Default is true.
%    'restrictBySd'      -Logical. If true, the parameter search restricts
%                         all params to within three standard deviations of
%                         their means. Default is true.
%    testingParamToVary2 -Character vector of an additional parameter to
%                         vary across groups of observers. Current options
%                         include 'noise' and'nMatches'. Default is [].
%    testingValsToVary2  -Vector which includes the different values of the
%                         second variable parameter to test. Default is [].
%    standardizeNMatches -Logical to be used with the 'varyTestWlIncr'
%                         method. When true, keeps the number of matches
%                         constant by adding multiple matches per test
%                         wavelength when fewer wavelengths are used (note
%                         that this is not exact, due to rounding).
%                         Default is false.
%    restrictErrWls      -Logical to be used with the 'varyTestWlIncr'
%                         method. When true, calculates average cone and
%                         match error using only wavelengths from the
%                         session with the least number of test
%                         wavelengths. Default is true.
%    freezeNoise         -Logical. Whe true, uses the same sampled
%                         observers for all match sets. Default is true.

% History:
%   07/23/20  dce       Wrote it.
%   07/27/20  dce       Added option to vary two parameters
%   07/28/20  dce       Cleaned up, added option to standardize number of
%                       wavelengths, added timing information
%   07/29/20  dce       Added noise freezing option
%   07/31/20  dce       Restructured the code
%   08/05/20  dce       Modified to calculate opponent contrast params

% Example:
%   sampleRayleighMatchSeries('test100',20,670,560,570:5:640,...
%   [0 0 1 1 0 1 1], 0.02,'noise', [0 0.01 0.02 0.04])
 
% Pause if there are any errors - helps with debugging
dbstop if error;

% Input parsing
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('baseConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('opponentParams',[0.8078 4.1146 1.2592 0.0200],@(x)(isvector(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',3201,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('noiseScaleFactor',0,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('testingParamToVary2',[],@(x)(ischar(x)));
p.addParameter('testingValsToVary2',[],@(x)(isnumeric(x)));
p.addParameter('standardizeNMatches',false,@(x)(islogical(x)));
p.addParameter('restrictErrWls',true,@(x)(islogical(x)));
p.addParameter('freezeNoise',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Input check
if isempty(p.Results.testingParamToVary2) ~= isempty(p.Results.testingValsToVary2)
    error('If a second variable parameter is specified, its values must also be specified');
end

% Set up directory for saving results
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearchSeries',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Make data-storing arrays. These are one-dimensional if one parameter is
% being varied, and two-dimensional otherwise.
if isempty(p.Results.testingParamToVary2)
    coneErrFC = zeros(length(testingValsToVary),1);      % Forced-choice
    coneErrStdFC = zeros(length(testingValsToVary),1); 
    matchErrFC = zeros(length(testingValsToVary),1);
    matchErrStdFC = zeros(length(testingValsToVary),1);
    matchErrSampledFC = zeros(length(testingValsToVary),1);
    coneErrAdjust = zeros(length(testingValsToVary),1);   % Adjustment
    coneErrStdAdjust = zeros(length(testingValsToVary),1);
    matchErrAdjust = zeros(length(testingValsToVary),1);
    matchErrStdAdjust = zeros(length(testingValsToVary),1);
    matchErrSampledAdjust = zeros(length(testingValsToVary),1);
else
    coneErrFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    coneErrAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    coneErrStdFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrStdFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    coneErrStdAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrStdAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrSampledFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrSampledAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
end

% Setup for collecting filenames
fileNames = cell(4*numel(matchErrFC),1);
counter = 1;

% Freeze noise if desired. This means that the same sampled observers will
% be used in all trials
if p.Results.freezeNoise
    sampledObservers = sampleRayleighObservers(nObservers,...
        p.Results.baseConeParams,coneParamsToVary);
else
    sampledObservers = [];
end

% When varying wl increment, calculate match errors using only the
% wavelengths from the coarsest match, if desired
if strcmp(testingParamToVary,'testWlIncr')
    if p.Results.restrictErrWls 
        errWls = test(1):max(testingValsToVary):test(end);
    else
        errWls = [];
    end
end
        
% Define the variable parameter names before starting loop 
if strcmp(testingParamToVary,'testWlIncr')...
        && strcmp(p.Results.testingParamToVary2,'nMatches')
    paramName = 'Number of Test Wavelengths';
    paramName2 = 'Number of Matches Per Test Light';
elseif strcmp(testingParamToVary,'testWlIncr')...
        && strcmp(p.Results.testingParamToVary2,'noise')
    paramName = 'Number of Test Wavelengths';
    paramName2 = 'Observer Noise';
elseif strcmp(testingParamToVary,'nMatches')...
        && strcmp(p.Results.testingParamToVary2,'noise')
    paramName = 'Number of Matches Per Test Light';
    paramName2 = 'Observer Noise';
elseif strcmp(testingParamToVary,'testWlIncr')
    if p.Results.standardizeNMatches
        nWls = ceil(length(test(1):test(end))*(1./testingValsToVary));
        matchesPerTrial = max(nWls);
        nMatches = floor(matchesPerTrial./nWls);
        paramName = 'Number of Test Wavelengths (standardized across trials)';
    else
        nMatches = p.Results.nObserverMatches*ones(1,length(testingValsToVary));
        paramName = 'Number of Test Wavelengths';
    end
elseif strcmp(testingParamToVary,'noise')
    paramName = 'Observer Noise Standard Deviation';
elseif strcmp(testingParamToVary,'nMatches')
    paramName = 'Number of Matches Per Test Light';
else
    error('Unrecognized parameter combination');
end
%% Run the simulation in a loop
% Vary wavelength increment and number of matches per test wavelength
for i = 1:length(testingValsToVary)
    if ~isempty(p.Results.testingParamToVary2) % Varying 2 params
        for j = 1:length(p.Results.testingValsToVary2)
            trialID = [subjID '_' num2str(testingValsToVary(i)) '_'...
                num2str(p.Results.testingValsToVary2(j))];
            
            % Vary wavelength increment and nMatches
            if strcmp(testingParamToVary,'testWlIncr')...
                    && strcmp(p.Results.testingParamToVary2,'nMatches')
                % Define parameters
                testSpds = test(1):testingValsToVary(i):test(end);
                nObserverMatches = p.Results.testingValsToVary2(j);
                plotTitle = ['nMatches = ' num2str(nObserverMatches)...
                    ', nWavelengths = ' num2str(length(testSpds))];
                % Run the simulation in a loop
                [coneErrAdjust(i,j),matchErrAdjust(i,j),coneErrStdAdjust(i,j),...
                    matchErrStdAdjust(i,j),matchErrSampledAdjust(i,j)] = ...
                    sampleRayleighMatch([trialID '_threshold'],nObservers,...
                    p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                    p1,p2,testSpds,'threshold','nObserverMatches',...
                    nObserverMatches,'nReversals',p.Results.nReversals,...
                    'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
                    'nBelowThreshold',p.Results.nBelowThreshold,'age',...
                    p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                    p.Results.p2Scale,'testScale',p.Results.testScale,...
                    'LMEqualOD',p.Results.LMEqualOD,'dlens0',p.Results.dlens0,...
                    'restrictBySd',p.Results.restrictBySd,'dmac0',...
                    p.Results.dmac0,'plotTitle',[plotTitle ', Threshold'],...
                    'adjustmentLength',p.Results.adjustmentLength,...
                    'sampledObservers',sampledObservers,'noiseScaleFactor',...
                    p.Results.noiseScaleFactor,'fieldSize',...
                    p.Results.fieldSize,'errWls',errWls);
                [coneErrFC(i,j),matchErrFC(i,j),coneErrStdFC(i,j),...
                    matchErrStdFC(i,j),matchErrSampledFC(i,j)] = ...
                    sampleRayleighMatch([trialID '_FC'],nObservers,...
                    p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                    p1,p2,testSpds,'forcedChoice','nObserverMatches',...
                    nObserverMatches,'nReversals',p.Results.nReversals,...
                    'age',p.Results.age,'p1Scale',p.Results.p1Scale,...
                    'p2Scale',p.Results.p2Scale,'testScale',....
                    p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                    'dlens0',p.Results.dlens0,'restrictBySd',...
                    p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                    'adjustmentLength',p.Results.adjustmentLength,...
                    'sampledObservers',sampledObservers,'plotTitle',...
                    [plotTitle ', FC'],'noiseScaleFactor',...
                    p.Results.noiseScaleFactor,'fieldSize',...
                    p.Results.fieldSize,'errWls',errWls);
                
                % Vary wavelength increment and noise
            elseif strcmp(testingParamToVary,'testWlIncr')...
                    && strcmp(p.Results.testingParamToVary2,'noise')
                % Define the parameters
                testSpds = test(1):testingValsToVary(i):test(end);
                plotTitle = ['Noise = ' num2str(p.Results.testingValsToVary2(j)...
                    *opponentNoise) ', nWavelengths = ' num2str(length(testSpds))];
                % Run the simulation with the two methods
                [coneErrAdjust(i,j),matchErrAdjust(i,j),coneErrStdAdjust(i,j),...
                    matchErrStdAdjust(i,j),matchErrSampledAdjust(i,j)] = ...
                    sampleRayleighMatch([trialID '_threshold'],nObservers,...
                    p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                    p1,p2,testSpds,'threshold','nObserverMatches',...
                    p.Results.nObserverMatches,'fieldSize',p.Results.fieldSize,...
                    'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
                    p.Results.thresholdScaleFactor,'nBelowThreshold',...
                    p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
                    p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
                    p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                    'dlens0',p.Results.dlens0,'restrictBySd',...
                    p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                    'adjustmentLength',p.Results.adjustmentLength,...
                    'sampledObservers',sampledObservers,'plotTitle',...
                    [plotTitle ', Threshold'],'noiseScaleFactor',...
                    p.Results.testingValsToVary2(j),'errWls',errWls);
                [coneErrFC(i,j),matchErrFC(i,j),coneErrStdFC(i,j),...
                    matchErrStdFC(i,j),matchErrSampledFC(i,j)] = ...
                    sampleRayleighMatch([trialID '_FC'],nObservers,...
                    p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                    p1,p2,testSpds,'forcedChoice','nObserverMatches',...
                    p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
                    'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                    p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
                    p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
                    p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                    'adjustmentLength',p.Results.adjustmentLength,...
                    'sampledObservers',sampledObservers,'plotTitle',...
                    [plotTitle ', FC'],'errWls',errWls,'noiseScaleFactor',...
                    p.Results.testingValsToVary2(j),'fieldSize',p.Results.fieldSize);
                
                % Vary number of matches and noise
            elseif strcmp(testingParamToVary,'nMatches')...
                    && strcmp(p.Results.testingParamToVary2,'noise')
                % Define the parameters
                nObserverMatches = testingValsToVary(i);
                plotTitle = ['Noise = ' num2str(p.Results.testingValsToVary2(j)...
                    *opponentNoise) ', nMatches = ' num2str(nObserverMatches)];
                % Run the simulation with the two methods
                [coneErrAdjust(i,j),matchErrAdjust(i,j),coneErrStdAdjust(i,j),...
                    matchErrStdAdjust(i,j),matchErrSampledAdjust(i,j)] = ...
                    sampleRayleighMatch([trialID '_threshold'],nObservers,...
                    p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                    p1,p2,test,'threshold','nObserverMatches',nObserverMatches,...
                    'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
                    p.Results.thresholdScaleFactor,'nBelowThreshold',...
                    p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
                    p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
                    p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                    'dlens0',p.Results.dlens0,'restrictBySd',...
                    p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                    'adjustmentLength',p.Results.adjustmentLength,...
                    'sampledObservers',sampledObservers,'plotTitle',...
                    [plotTitle ', Threshold'],'noiseScaleFactor',...
                    p.Results.testingValsToVary2(j),'fieldSize',p.Results.fieldSize);
                [coneErrFC(i,j),matchErrFC(i,j),coneErrStdFC(i,j),...
                    matchErrStdFC(i,j),matchErrSampledFC(i,j)] = ...
                    sampleRayleighMatch([trialID '_FC'],nObservers,...
                    p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                    p1,p2,test,'forcedChoice','nObserverMatches',...
                    nObserverMatches,'nReversals',p.Results.nReversals,...
                    'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                    p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
                    p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
                    p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                    'adjustmentLength',p.Results.adjustmentLength,...
                    'sampledObservers',sampledObservers,'plotTitle',...
                    [plotTitle ', FC'],'noiseScaleFactor',...
                    p.Results.testingValsToVary2(j),'fieldSize',p.Results.fieldSize);
            end
            % Save filenames
            fileNames(counter) = cellstr(fullfile(getpref('ForcedChoiceCM',...
                'rayleighDataDir'),'paramsSearch',[trialID '_threshold'],...
                [trialID '_threshold_paramPlots.pdf']));
            fileNames(counter+1) = cellstr(fullfile(getpref('ForcedChoiceCM',...
                'rayleighDataDir'),'paramsSearch',[trialID '_threshold']...
                ,[trialID '_threshold_errPlots.pdf']));
            fileNames(counter+2) = cellstr(fullfile(getpref('ForcedChoiceCM',...
                'rayleighDataDir'),'paramsSearch',[trialID '_FC'],...
                [trialID '_FC_paramPlots.pdf']));
            fileNames(counter+3) = cellstr(fullfile(getpref('ForcedChoiceCM',...
                'rayleighDataDir'),'paramsSearch',[trialID '_FC']...
                ,[trialID '_FC_errPlots.pdf']));
            counter = counter+4;
            save(fullfile(outputDir,'testData.mat'),'p','coneErrFC',...
                'matchErrFC','coneErrAdjust','matchErrAdjust',...
                'coneErrStdFC','matchErrStdFC','matchErrSampledFC',...
                'coneErrStdAdjust','matchErrStdAdjust',...
                'matchErrSampledAdjust','subjID','nObservers','p1','p2',...
                'test','coneParamsToVary','testingParamToVary','outputDir',...
                'testingValsToVary','paramName','paramName2','fileNames');
        end
    else  % Varying one param
        trialID = [subjID '_' num2str(testingValsToVary(i))];
        % Vary test wavelength increment
        if strcmp(testingParamToVary,'testWlIncr')
            % Define the range of test wavelengths
            testSpds = test(1):testingValsToVary(i):test(end);
            plotTitle = ['nWavelengths = ' num2str(length(testSpds))];
            % Run the simulation with the two methods
            [coneErrAdjust(i),matchErrAdjust(i),coneErrStdAdjust(i),...
                matchErrStdAdjust(i),matchErrSampledAdjust(i)] = ...
                sampleRayleighMatch([trialID '_threshold'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                p1,p2,testSpds,'threshold','nObserverMatches',nMatches(i),...
                'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
                p.Results.thresholdScaleFactor,'nBelowThreshold',...
                p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
                p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
                p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                'dlens0',p.Results.dlens0,'restrictBySd',...
                p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers,'plotTitle',...
                [plotTitle ', Threshold'],'errWls',errWls,'noiseScaleFactor',...
                p.Results.noiseScaleFactor,'fieldSize',p.Results.fieldSize);
            [coneErrFC(i),matchErrFC(i),coneErrStdFC(i),...
                matchErrStdFC(i),matchErrSampledFC(i)] = ...
                sampleRayleighMatch([trialID '_FC'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                p1,p2,testSpds,'forcedChoice','nObserverMatches',nMatches(i),...
                'nReversals',p.Results.nReversals,'age',p.Results.age,...
                'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,...
                'testScale',p.Results.testScale,'LMEqualOD',...
                p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
                p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers,'plotTitle',...
                [plotTitle ', FC'],'errWls',errWls,'noiseScaleFactor',...
                p.Results.noiseScaleFactor,'fieldSize',p.Results.fieldSize);
            
            % Vary noise
        elseif strcmp(testingParamToVary,'noise')
            plotTitle = ['Noise = ' num2str(testingValsToVary(i))];
            % Sample, test, and recover parameters with the two methods
            [coneErrAdjust(i),matchErrAdjust(i),coneErrStdAdjust(i),...
                matchErrStdAdjust(i),matchErrSampledAdjust(i)] = ...
                sampleRayleighMatch([trialID '_threshold'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,p1,p2,test,...
                'threshold','nObserverMatches',p.Results.nObserverMatches,...
                'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
                p.Results.thresholdScaleFactor,'nBelowThreshold',...
                p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
                p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
                p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                'dlens0',p.Results.dlens0,'restrictBySd',p.Results.restrictBySd,...
                'dmac0',p.Results.dmac0,'adjustmentLength',....
                p.Results.adjustmentLength,'sampledObservers',sampledObservers,...
                'plotTitle',[plotTitle ', Threshold'],'noiseScaleFactor',...
                testingValsToVary(i),'fieldSize',p.Results.fieldSize);
            [coneErrFC(i),matchErrFC(i),coneErrStdFC(i),matchErrStdFC(i),...
                matchErrSampledFC(i)] = ...
                sampleRayleighMatch([trialID '_FC'],nObservers,...
                p.Results.baseConeParams,...
                coneParamsToVary,p.Results.opponentParams,p1,p2,test,'forcedChoice',...
                'nObserverMatches',p.Results.nObserverMatches,'nReversals',...
                p.Results.nReversals,'age',p.Results.age,'p1Scale',...
                p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
                p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                'dlens0',p.Results.dlens0,'restrictBySd',...
                p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers, 'plotTitle',...
                [plotTitle ', FC'],'noiseScaleFactor',...
                testingValsToVary(i),'fieldSize',p.Results.fieldSize);
            
            % Vary number of matches per test light
        elseif strcmp(testingParamToVary,'nMatches')
            plotTitle = ['nMatches = ' num2str(testingValsToVary(i))];
            [coneErrAdjust(i),matchErrAdjust(i),coneErrStdAdjust(i),...
                matchErrStdAdjust(i),matchErrSampledAdjust(i)] = ...
                sampleRayleighMatch([trialID '_threshold'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,...
                p1,p2,test,'threshold','nObserverMatches',...
                testingValsToVary(i),'nReversals',p.Results.nReversals,...
                'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
                'nBelowThreshold',p.Results.nBelowThreshold,'age',...
                p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                p.Results.p2Scale,'testScale',p.Results.testScale,...
                'LMEqualOD',p.Results.LMEqualOD,'dlens0',p.Results.dlens0,...
                'restrictBySd',p.Results.restrictBySd,'dmac0',...
                p.Results.dmac0,'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers,'plotTitle',...
                [plotTitle ', Threshold'],'noiseScaleFactor',...
                p.Results.noiseScaleFactor,'fieldSize',p.Results.fieldSize);
            [coneErrFC(i),matchErrFC(i),coneErrStdFC(i),...
                matchErrStdFC(i),matchErrSampledFC(i)] = ...
                sampleRayleighMatch([trialID '_FC'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p.Results.opponentParams,....
                p1,p2,test,'forcedChoice','nObserverMatches',...
                testingValsToVary(i),'nReversals',p.Results.nReversals,'age',...
                p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                p.Results.p2Scale,'testScale',p.Results.testScale,...
                'LMEqualOD',p.Results.LMEqualOD,'dlens0',p.Results.dlens0,...
                'restrictBySd',p.Results.restrictBySd,'dmac0',...
                p.Results.dmac0,'adjustmentLength',...
                p.Results.adjustmentLength,'sampledObservers',sampledObservers,...
                'plotTitle',[plotTitle ', FC'],'noiseScaleFactor',...
                p.Results.noiseScaleFactor,'fieldSize',p.Results.fieldSize);
        end
        % Save filenames
        fileNames(counter) = cellstr(fullfile(getpref('ForcedChoiceCM',...
            'rayleighDataDir'),'paramsSearch',[trialID '_threshold'],...
            [trialID '_threshold_paramPlots.pdf']));
        fileNames(counter+1) = cellstr(fullfile(getpref('ForcedChoiceCM',...
            'rayleighDataDir'),'paramsSearch',[trialID '_threshold']...
            ,[trialID '_threshold_errPlots.pdf']));
        fileNames(counter+2) = cellstr(fullfile(getpref('ForcedChoiceCM',...
            'rayleighDataDir'),'paramsSearch',[trialID '_FC'],...
            [trialID '_FC_paramPlots.pdf']));
        fileNames(counter+3) = cellstr(fullfile(getpref('ForcedChoiceCM',...
            'rayleighDataDir'),'paramsSearch',[trialID '_FC']...
            ,[trialID '_FC_errPlots.pdf']));
        counter = counter+4;
        save(fullfile(outputDir,'testData.mat'),'p','coneErrFC',...
            'matchErrFC','coneErrAdjust','matchErrAdjust','outputDir',...
            'coneErrStdFC','matchErrStdFC','matchErrSampledFC',...
            'coneErrStdAdjust','matchErrStdAdjust','matchErrSampledAdjust',...
            'subjID','nObservers','p1','p2','test','coneParamsToVary',...
            'testingParamToVary','testingValsToVary','paramName',...
            'fileNames');
    end
end

%% Plotting
% y limits for cone/match error plots
coneErrLim = 0.02;
matchErrLim = 0.2;

% Offsets for text labels
dx =  -0.001* [0,ones(1,length(testingValsToVary)-1)];
dyFC = 0.2;
dyThreshold = 0.3;

if isempty(p.Results.testingParamToVary2)% Make plots for a single variable
    theFig = figure();
    set(theFig,'Color',[1 1 1],'Position',[10 10 1400 800]);
    hold on;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 2, ...
        'colsNum', 1, ...
        'heightMargin',  0.1, ...
        'widthMargin',    0.1, ...
        'leftMargin',     0.1, ...
        'rightMargin',    0.1, ...
        'bottomMargin',   0.1, ...
        'topMargin',      0.1);
    
    % Subplot 1 - cone spectral sensitivity error
    subplot('Position', subplotPosVectors(1,1).v);
    hold on;
    % If wavelength increment is being varied, the x value is the number of
    % wavelengths, not the increment.
    if strcmp(testingParamToVary,'testWlIncr')
        xVals = ceil(length(570:610)*(1./testingValsToVary));
        %xVals = ceil(length(test(1):test(end))*(1./testingValsToVary));
    else
        xVals = testingValsToVary;
    end
    % Plot results
    plot(xVals,coneErrFC,'r-o','LineWidth',2.5);
    plot(xVals,coneErrStdFC,'g-o','LineWidth',2.5);
    if ~isempty(coneErrAdjust)
        plot(xVals,coneErrAdjust,'b-o','LineWidth',2.5);
        labels2 = cellstr(num2str(coneErrAdjust,3));
        text(xVals+dx,coneErrFC+dyThreshold*coneErrLim,labels2,'Color','b');
        legend('Forced Choice - Recovered Cones',...
            'Forced Choice - Standard Cones','Adjustment - Recovered Cones');
    else
        legend('Forced Choice','Forced Choice - Standard Cones');
    end
    labels1 = cellstr(num2str(coneErrFC,3));
    text(xVals+dx,coneErrFC+dyFC*coneErrLim,labels1,'Color','r');
    labels3 = cellstr(num2str(coneErrStdFC,3));
    text(xVals+dx,coneErrStdFC+dyFC*coneErrLim*0.5,labels3,'Color','g');
    
    title('Average Cone Spectral Sensitivity Error');
    xlabel(paramName);
    ylabel('Average RMS Error');
    ylim([0 coneErrLim]);
    text(xVals(end)/10,0.0175,['n = ' num2str(nObservers)]);
    
    % Subplot 2 - Match error
    subplot('Position', subplotPosVectors(2,1).v);
    hold on;
    plot(xVals,matchErrFC,'r-o','LineWidth',2.5);
    plot(xVals,matchErrStdFC,'g-o','LineWidth',2.5);
    plot(xVals,matchErrSampledFC,'y-o','LineWidth',2.5);
    if ~strcmp(testingParamToVary,'adjustmentLength')
        plot(xVals,matchErrAdjust,'b-o','LineWidth',2.5);
        labels2 = cellstr(num2str(matchErrAdjust,3));
        text(xVals+dx,matchErrFC+dyThreshold*matchErrLim,labels2,'Color','b')
        legend('Forced Choice','Forced Choice - Standard Cones',...
            'Forced Choice - Sampled Cones','Adjustment');
    else
        legend('Forced Choice','Forced Choice - Standard Cones',...
            'Forced Choice - Sampled Cones');
    end
    labels1 = cellstr(num2str(matchErrFC,3));
    text(xVals+dx,matchErrFC+dyFC*matchErrLim,labels1,'Color','r');
    labels3 = cellstr(num2str(matchErrStdFC,3));
    text(xVals+dx,matchErrStdFC+dyFC*matchErrLim*0.5,labels3,'Color','g');
    labels4 = cellstr(num2str(matchErrSampledFC,3));
    text(xVals+dx,matchErrFC+0.5*dyFC*matchErrLim,labels4,'Color','y');
    
    title('Average Match Error');
    xlabel(paramName);
    ylabel('Error');
    ylim([0 matchErrLim]);
    text(xVals(end)/10,0.045,['n = ' num2str(nObservers)]);
    sgtitle(['Vary ' paramName]);
    
else            % Make plots for variation of two parameters
    theFig = figure();
    set(theFig,'Color',[1 1 1],'Position',[10 10 1400 800]);
    hold on;
    subplotPosVectors = NicePlot.getSubPlotPosVectors(...
        'rowsNum', 2, ...
        'colsNum', 2, ...
        'heightMargin',  0.1, ...
        'widthMargin',    0.1, ...
        'leftMargin',     0.1, ...
        'rightMargin',    0.1, ...
        'bottomMargin',   0.1, ...
        'topMargin',      0.1);
    
    % If wavelength increment is being varied, the y value is the number of
    % wavelengths, not the increment.
    if strcmp(testingParamToVary,'testWlIncr')
        yVals = ceil(length(test(1):test(end))*(1./testingValsToVary));
    else
        yVals = testingValsToVary;
    end
    xVals = p.Results.testingValsToVary2;
    
    % Subplot 1 - cone spectral sensitivity error threshold
    subplot('Position', subplotPosVectors(1,1).v);
    theData = coneErrAdjust;
    heatmap(xVals,yVals,theData,'ColorLimits',[0 coneErrLim]);
    title('Cone Error - Threshold');
    ylabel(paramName);
    xlabel(paramName2);
    
    % Subplot 2 - Match error threshold
    subplot('Position', subplotPosVectors(2,1).v);
    theData = matchErrAdjust;
    heatmap(xVals,yVals,theData,'ColorLimits',[0 matchErrLim]);
    title('Match Error - Threshold');
    ylabel(paramName);
    xlabel(paramName2);
    
    % Subplot 3 - cone spectral sensitivity error forced choice
    subplot('Position', subplotPosVectors(1,2).v);
    theData = coneErrFC;
    heatmap(xVals,yVals,theData,'ColorLimits',[0 coneErrLim]);
    title('Cone Error - Forced Choice');
    ylabel(paramName);
    xlabel(paramName2);
    
    % Subplot 4 - Match error threshold
    subplot('Position', subplotPosVectors(2,2).v);
    theData = matchErrFC;
    heatmap(xVals,yVals,theData, 'ColorLimits',[0 matchErrLim]);
    title('Match Error - Forced Choice');
    ylabel(paramName);
    xlabel(paramName2);
    sgtitle(['Vary ' paramName ' and ' paramName2]);
end

%% Save plots
NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_plots']),theFig,300);

% Append the individual sampling plots to create one document
fileNames = string(fileNames);
append_pdfs(fullfile(outputDir,[subjID '_plots.pdf']),fileNames)
end