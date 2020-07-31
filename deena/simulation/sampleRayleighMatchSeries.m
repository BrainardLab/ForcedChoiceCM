function sampleRayleighMatchSeries(subjID,nObservers,p1,p2,test,...
    coneParamsToVary,testingParamToVary,testingValsToVary,varargin)
% Runs Rayleigh matching simulations with different groups of observers
% Syntax:
%   sampleRayleighMatchSeries(subjID,nObservers,p1,p2,test,coneParamsToVary,
%   testingParamToVary,testingValsToVary,varargin)
%
% Description:
%    This program tests how varying different experimental settings
%    affects the ability to recover observer parameters using
%    sampleRayleighMatch. The experimenter specifies a testing parameter to
%    systematically vary (currently, options include observer noise,
%    increment of test wavelengths, simulation adjustment length, and
%    number of matches per trial). For each value of the parameter, the
%    program runs sampleRayleighMatch, which samples observers, simulates
%    Rayleigh matching, and attempts to recover cone parameters. Both the
%    forced-choice and threshold versions are tested for each variable
%    parameter value (with the exception of simulation adjustment length,
%    which is only run for forced choice).
%
%    The program makes and saves plots of two types of error - average cone
%    spectral sensitivity error and average match error - allowing
%    comparison across different values of the variable parameter.
%
%    There are also options to vary test wavelength increment along with an
%    additional parameter (current options include noise or number of
%    matches per test wavelength). In this case, comparative "heatmaps" of
%    the two parameters are produced instead of the single parameter plots.
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
%                        'noise','nMatches', 'adjustmentLength', and
%                        'testWlIncr'.
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
%    'baseConeParams'    -Nine-element numeric vector of individual
%                         difference parameters used as a starting point.
%                         Default is zeros(1,9)
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
%    freezeNoise         -Logical. Whe true, uses the same sampled
%                         parameters for all match sets. Default is true.

% History:
%   07/23/20  dce       Wrote it.
%   07/27/20  dce       Added option to vary two parameters
%   07/28/20  dce       Cleaned up, added option to standardize number of
%                       wavelengths, added timing information
%   07/29/20  dce       Added noise freezing option
%   07/31/20  dce       Changed figure naming

% Example:
%   sampleRayleighMatchSeries('test100',20,670,560,570:5:640,...
%   [0 0 1 1 0 1 1 0], 'noise', [0 0.01 0.02 0.04])
dbstop if error; % Pause if there are any errors - helps with debugging
% Input parsing
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('baseConeParams',zeros(1,9),@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('adjustmentLength',3201,@(x)(isnumeric(x)));
p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',true,@(x)(islogical(x)));
p.addParameter('dmac0',true,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.addParameter('testingParamToVary2',[],@(x)(ischar(x)));
p.addParameter('testingValsToVary2',[],@(x)(isnumeric(x)));
p.addParameter('standardizeNMatches',false,@(x)(islogical(x)));
p.addParameter('freezeNoise',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Input checks
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
    matchErrFC = zeros(length(testingValsToVary),1);
    coneErrAdjust = zeros(length(testingValsToVary),1);  % Adjustment
    matchErrAdjust = zeros(length(testingValsToVary),1);
else
    coneErrFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrFC = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    coneErrAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
    matchErrAdjust = zeros(length(testingValsToVary),...
        length(p.Results.testingValsToVary2));
end
% Get ready to collect filenames
[row,col] = size(matchErrFC);
fileNames = cell(4*row*col,1);
counter = 1;

% Freeze noise if desired
if p.Results.freezeNoise
    sampledObservers = sampleRayleighObservers(nObservers,...
        p.Results.baseConeParams,coneParamsToVary);
else
    sampledObservers = [];
end

%% Define the variable parameter(s), then run the simulation in a loop
% Vary wavelength increment and number of matches per test wavelength
if strcmp(testingParamToVary,'testWlIncr')...
        && strcmp(p.Results.testingParamToVary2,'nMatches')
    paramName = 'Number of Test Wavelengths';
    paramName2 = 'Number of Matches Per Test Light';
    for i = 1:length(testingValsToVary)
        for j = 1:length(p.Results.testingValsToVary2)
            % Define the range of test wavelengths
            testSpds = test(1):testingValsToVary(i):test(end);
            nObserverMatches = p.Results.testingValsToVary2(j);
            trialID = [subjID '_' num2str(testingValsToVary(i)) '_'...
                num2str(p.Results.testingValsToVary2(j))];
            plotTitle = ['nMatches = ' num2str(nObserverMatches)...
                ', nWavelengths = ' num2str(length(testSpds))];
            
            % Run the simulation with the two methods
            [coneErrAdjust(i,j),matchErrAdjust(i,j)] = ...
                sampleRayleighMatch([trialID '_threshold'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p1,p2,testSpds,...
                'threshold','nObserverMatches',nObserverMatches,...
                'nReversals',p.Results.nReversals,...
                'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
                'nBelowThreshold',p.Results.nBelowThreshold,'age',...
                p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                p.Results.p2Scale,'testScale',p.Results.testScale,...
                'LMEqualOD',p.Results.LMEqualOD,'dlens0',p.Results.dlens0,...
                'restrictBySd',p.Results.restrictBySd,'dmac0',...
                p.Results.dmac0,'plotTitle',[plotTitle ' Threshold'],...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers);
            [coneErrAdjust(i,j),matchErrAdjust(i,j)] = ...
                sampleRayleighMatch([trialID '_FC'],nObservers,...
                p.Results.baseConeParams,coneParamsToVary,p1,p2,testSpds,...
                'forcedChoice','nObserverMatches',nObserverMatches,...
                'nReversals',p.Results.nReversals,'age',p.Results.age,...
                'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,...
                'testScale',p.Results.testScale,'LMEqualOD',...
                p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
                p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers,'plotTitle',...
                [plotTitle ' FC']);
            
            % Collect plot titles 
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
                'matchErrFC','coneErrAdjust','matchErrAdjust','subjID',...
                'nObservers','p1','p2','test','coneParamsToVary',...
                'testingParamToVary','testingValsToVary','paramName',...
                'paramName2','fileNames');
        end
    end
    
    % Vary wavelength increment and noise
elseif strcmp(testingParamToVary,'testWlIncr')...
        && strcmp(p.Results.testingParamToVary2, 'noise')
    paramName = 'Number of Test Wavelengths';
    paramName2 = 'Observer Noise';
    for i = 1:length(testingValsToVary)
        for j = 1:length(p.Results.testingValsToVary2)
            % Define the range of test wavelengths
            testSpds = test(1):testingValsToVary(i):test(end);
            trialConeParams = [p.Results.baseConeParams(1:8),...
                p.Results.testingValsToVary2(j)];
            if p.Results.freezeNoise
                sampledObservers(:,9) = p.Results.testingValsToVary2(j);
            end 
            trialID = [subjID '_' num2str(testingValsToVary(i)) '_'...
                num2str(p.Results.testingValsToVary2(j))];
            plotTitle = ['noise = ' num2str(p.Results.testingValsToVary2(j))...
                ', nWavelengths = ' num2str(length(testSpds))];
            
            % Run the simulation with the two methods
            [coneErrAdjust(i,j),matchErrAdjust(i,j)] = ...
                sampleRayleighMatch([trialID '_threshold'],nObservers,...
                trialConeParams,coneParamsToVary,p1,p2,testSpds,...
                'threshold','nObserverMatches',p.Results.nObserverMatches,...
                'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
                p.Results.thresholdScaleFactor,'nBelowThreshold',...
                p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
                p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
                p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
                'dlens0',p.Results.dlens0,'restrictBySd',...
                p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers,'plotTitle',...
                [plotTitle ' Threshold']);
            [coneErrFC(i,j),matchErrFC(i,j)] = ...
                sampleRayleighMatch([trialID '_FC'],nObservers,...
                trialConeParams,coneParamsToVary,p1,p2,testSpds,...
                'forcedChoice','nObserverMatches',...
                p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
                'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
                p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
                p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
                p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
                'adjustmentLength',p.Results.adjustmentLength,...
                'sampledObservers',sampledObservers,'plotTitle',...
                [plotTitle ' FC']);
            
            % Save plot filenames
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
                'rayleighDataDir'),'paramsSearch',[trialID 'FC']...
                ,[trialID '_FC_errPlots.pdf']));
            counter = counter+4;
            save(fullfile(outputDir,'testData.mat'),'p','coneErrFC',...
                'matchErrFC','coneErrAdjust','matchErrAdjust','subjID',...
                'nObservers','p1','p2','test','coneParamsToVary',...
                'testingParamToVary','testingValsToVary','paramName',...
                'paramName2','fileNames');
        end
    end
    
    % Vary wavelength increment
elseif strcmp(testingParamToVary,'testWlIncr')
    if p.Results.standardizeNMatches
        nWls = ceil(length(test(1):test(end))*(1./testingValsToVary));
        matchesPerTrial = max(nWls);
        nMatches = floor(matchesPerTrial./nWls);
        paramName = 'Number of Test Wavelengths (standardized across trials)';
    else
        nMatches = p.Results.nMatches*ones(1,length(testingValsToVary));
        paramName = 'Number of Test Wavelengths';
    end
    
    for i = 1:length(testingValsToVary)
        % Define the range of test wavelengths
        testSpds = test(1):testingValsToVary(i):test(end);
        trialID = [subjID '_' num2str(testingValsToVary(i))];
        plotTitle = ['nWavelengths = ' num2str(length(testSpds))];
        
        % Run the simulation with the two methods
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
            sampleRayleighMatch([trialID '_threshold'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,testSpds,...
            'threshold','nObserverMatches',nMatches(i),...
            'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
            p.Results.thresholdScaleFactor,'nBelowThreshold',...
            p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
            p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
            p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
            'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength,...
            'sampledObservers',sampledObservers,'plotTitle',...
            [plotTitle ' Threshold']);
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,testSpds,...
            'forcedChoice','nObserverMatches',nMatches(i),...
            'nReversals',p.Results.nReversals,'age',p.Results.age,...
            'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,...
            'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength,...
            'sampledObservers',sampledObservers,'plotTitle',...
            [plotTitle ' FC']);
        
        % Save plot filenames 
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
            'rayleighDataDir'),'paramsSearch',[trialID 'FC']...
            ,[trialID '_FC_errPlots.pdf']));
        counter = counter+4;     
        save(fullfile(outputDir,'testData.mat'),'p','coneErrFC',...
            'matchErrFC','coneErrAdjust','matchErrAdjust','subjID',...
            'nObservers','p1','p2','test','coneParamsToVary',...
            'testingParamToVary','testingValsToVary','paramName',...
            'fileNames');
    end
    
    % Vary noise
elseif strcmp(testingParamToVary,'noise')
    paramName = 'Observer Noise Standard Deviation';
    for i = 1:length(testingValsToVary)
        % Make base params with the specified noise level
        trialConeParams = [p.Results.baseConeParams(1:8),...
            testingValsToVary(i)];
        if p.Results.freezeNoise
            sampledObservers(:,9) = testingValsToVary(i);
        end
        trialID = [subjID '_' num2str(testingValsToVary(i))];
        plotTitle = ['noise = ' num2str(testingValsToVary(i))];  
        % Sample, test, and recover parameters with the two methods
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
            sampleRayleighMatch([trialID '_threshold'],nObservers,...
            trialConeParams,coneParamsToVary,p1,p2,test,...
            'threshold','nObserverMatches',p.Results.nObserverMatches,...
            'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
            p.Results.thresholdScaleFactor,'nBelowThreshold',...
            p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
            p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
            p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
            'dlens0',p.Results.dlens0,'restrictBySd',p.Results.restrictBySd,...
            'dmac0',p.Results.dmac0,'adjustmentLength',....
            p.Results.adjustmentLength,'sampledObservers',sampledObservers,...
            'plotTitle',[plotTitle ' Threshold']);          
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,trialConeParams,...
            coneParamsToVary,p1,p2,test,'forcedChoice','nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength,...
            'sampledObservers',sampledObservers, 'plotTitle',...
            [plotTitle ' FC']);
        
        % Save plot filenames 
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
            'matchErrFC','coneErrAdjust','matchErrAdjust','subjID',...
            'nObservers','p1','p2','test','coneParamsToVary',...
            'testingParamToVary','testingValsToVary','paramName',...
            'fileNames');
    end
    
    % Vary number of matches
elseif strcmp(testingParamToVary,'nMatches')
    paramName = 'Number of Matches Per Test Light';
    for i = 1:length(testingValsToVary)
        trialID = [subjID '_' num2str(testingValsToVary(i))];
        plotTitle = ['nMatches = ' num2str(testingValsToVary(i))];
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
            sampleRayleighMatch([trialID '_threshold'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,test,...
            'threshold','nObserverMatches',testingValsToVary(i),...
            'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
            p.Results.thresholdScaleFactor,'nBelowThreshold',...
            p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
            p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
            p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,'dlens0',...
            p.Results.dlens0,'restrictBySd',p.Results.restrictBySd,'dmac0',...
            p.Results.dmac0,'adjustmentLength',p.Results.adjustmentLength,...
            'sampledObservers',sampledObservers,'plotTitle',...
            [plotTitle ' Threshold']);
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,test,...
            'forcedChoice','nObserverMatches',testingValsToVary(i),...
            'nReversals',p.Results.nReversals,'age',p.Results.age,...
            'p1Scale',p.Results.p1Scale,'p2Scale',p.Results.p2Scale,...
            'testScale',p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
            'dlens0',p.Results.dlens0,'restrictBySd',p.Results.restrictBySd,...
            'dmac0',p.Results.dmac0,'adjustmentLength',...
            p.Results.adjustmentLength,'sampledObservers',sampledObservers,...
            'plotTitle',[plotTitle ' FC']);
        
        % Save plot filenames 
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
            'rayleighDataDir'),'paramsSearch',[trialID 'FC']...
            ,[trialID '_FC_errPlots.pdf'])); 
        save(fullfile(outputDir,'testData.mat'),'p','coneErrFC',...
            'matchErrFC','coneErrAdjust','matchErrAdjust','subjID',...
            'nObservers','p1','p2','test','coneParamsToVary',...
            'testingParamToVary','testingValsToVary','paramName',...
            'fileNames');
    end
    
    % Vary simulation adjustment length. Note that this only runs the
    % forced choice simulation, as threshold often cannot be reached for
    % coarser arrays
elseif strcmp(testingParamToVary,'adjustmentLength')
    paramName = 'Number of Possible Light Adjustments';
    fileNames = cell(2*row*col,1);
    for i = 1:length(testingValsToVary)
        trialID = [subjID '_' num2str(testingValsToVary(i))];
        plotTitle = ['adjustmentLength = ' num2str(testingValsToVary(i))];
        
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,test,...
            'forcedChoice','nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',testingValsToVary(i),...
            'sampledObservers',sampledObservers,'plotTitle',[plotTitle ' FC']);
        
        % Save plot filenames 
        fileNames(counter) = cellstr(fullfile(getpref('ForcedChoiceCM',...
            'rayleighDataDir'),'paramsSearch',[trialID '_threshold'],...
            [trialID '_FC_paramPlots.pdf']));
        fileNames(counter+1) = cellstr(fullfile(getpref('ForcedChoiceCM',...
            'rayleighDataDir'),'paramsSearch',[trialID '_threshold']...
            ,[trialID '_FC_errPlots.pdf']));
        counter = counter+2;
        
        save(fullfile(outputDir,'testData.mat'),'p','coneErrFC',...
            'matchErrFC','coneErrAdjust','matchErrAdjust','subjID',...
            'nObservers','p1','p2','test','coneParamsToVary',...
            'testingParamToVary','testingValsToVary','paramName',...
            'fileNames');
    end
else
    error('Unrecognized parameter to vary');
end

%% Plotting
% y limits for cone/match error plots
coneErrLim = 0.01;
matchErrLim = 0.05;

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
        xVals = ceil(length(test(1):test(end))*(1./testingValsToVary));
    else
        xVals = testingValsToVary;
    end
    % Plot results
    plot(xVals,coneErrFC,'r-o','LineWidth',2.5);
    if ~isempty(coneErrAdjust)
        plot(xVals,coneErrAdjust,'b-o','LineWidth',2.5);
        labels2 = cellstr(num2str(coneErrAdjust,3));
        text(xVals+dx,coneErrFC+dyThreshold*coneErrLim,labels2,'Color','b');
        legend('Forced Choice','Adjustment');
    end
    labels1 = cellstr(num2str(coneErrFC,3));
    text(xVals+dx,coneErrFC+dyFC*coneErrLim,labels1,'Color','r');
    
    title('Average Cone Spectral Sensitivity Error');
    xlabel(paramName);
    ylabel('Average RMS Error');
    ylim([0 coneErrLim]);
    text(xVals(end)/10,0.0175,['n = ' num2str(nObservers)]);
    
    % Subplot 2 - Match error
    subplot('Position', subplotPosVectors(2,1).v);
    hold on;
    plot(xVals,matchErrFC,'r-o','LineWidth',2.5);
    if ~strcmp(testingParamToVary,'adjustmentLength')
        plot(xVals,matchErrAdjust,'b-o','LineWidth',2.5);
        labels2 = cellstr(num2str(matchErrAdjust,3));
        text(xVals+dx,matchErrFC+dyThreshold*matchErrLim,labels2,'Color','b')
        legend('Forced Choice','Adjustment');
    end
    labels1 = cellstr(num2str(matchErrFC,3));
    text(xVals+dx,matchErrFC+dyFC*matchErrLim,labels1,'Color','r');
    
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

% Save plots
NicePlot.exportFigToPDF(fullfile(outputDir,[subjID '_plots']),theFig,300);

% Append the individual sampling plots to create one document
fileNames = string(fileNames); % Convert files to string array
append_pdfs(fullfile(outputDir,[subjID '_plots.pdf']),fileNames)
end