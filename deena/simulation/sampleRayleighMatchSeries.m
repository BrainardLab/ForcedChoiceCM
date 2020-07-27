function sampleRayleighMatchSeries(subjID,nObservers,p1,p2,test,...
    coneParamsToVary,testingParamToVary,testingValsToVary,varargin)
% Runs Rayleigh matching simulations with different groups of observers

% Syntax:
%   sampleRayleighMatchSeries(subjID,nObservers,p1,p2,test,comeParamsToVary,
%   testingParamsToVary,testingValsToVary,varargin)
%
% Description:
%    This program tests how varying different experimental parameters
%    affects the ability to recover observer parameters using
%    sampleRayleighMatch. The experimenter specifies a testing parameter to
%    systematically vary (currently, options include observer noise,
%    increment of test wavelengths, simulation adjustment length, and 
%    number of matches per trial). For each value of the parameter, the 
%    program runs sampleRayleighMatch, which samples observers, simulates 
%    Rayleigh matching, and attempts to recover cone parameters. Both the 
%    forced-choice and threshold versions are tested for each variable 
%    parameter value. 
%
%    The program makes and saves plots of two types of error - average cone 
%    spectral sensitivity error and average match error - allowing 
%    comparison across different values of the variable parameter. 
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
 
% History:
%   07/23/20  dce       Wrote it.

% Example:
%   sampleRayleighMatchSeries('test100',20,670,560,570:5:640,...
%   [0 0 1 1 0 1 1 0], 'noise', [0 0.01 0.02 0.04])

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
p.parse(varargin{:});

% Make error-storing arrays
coneErrFC = zeros(length(testingValsToVary),1);      % Forced-choice
matchErrFC = zeros(length(testingValsToVary),1);
coneErrAdjust = zeros(length(testingValsToVary),1);  % Adjustment 
matchErrAdjust = zeros(length(testingValsToVary),1);

% Set up directory for saving results 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'paramsSearchSeries',subjID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Define the variable parameter, then run the simulation in a loop
% Vary noise
if strcmp(testingParamToVary,'noise')
    for i = 1:length(testingValsToVary)
        % Make base params with the specified noise level
        paramName = 'Observer Noise Standard Deviation';
        trialConeParams = [p.Results.baseConeParams(1:8),...
            testingValsToVary(i)]; 
        trialID = [subjID '_' num2str(i)];
        
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
        sampleRayleighMatch([trialID '_threshold'],nObservers,...
            trialConeParams,coneParamsToVary,p1,p2,test,...
            'threshold','plotResults',true,'nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
            'nBelowThreshold',p.Results.nBelowThreshold,'age',...
            p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength);
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,trialConeParams,...
            coneParamsToVary,p1,p2,test,'forcedChoice','plotResults',...
            true,'nObserverMatches',p.Results.nObserverMatches,...
            'nReversals',p.Results.nReversals,'age',p.Results.age,'p1Scale',...
            p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
            p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
            'dlens0',p.Results.dlens0,'restrictBySd',p.Results.restrictBySd,...
            'dmac0',p.Results.dmac0,'adjustmentLength',...
            p.Results.adjustmentLength);
    end
    
    % Vary wavelength increment
elseif strcmp(testingParamToVary,'testWlIncr')
    for i = 1:length(testingValsToVary)
        % Define the range of test wavelengths 
        paramName = 'Number of Test Wavelengths';
        testSpd = test(1):testingValsToVary(i):test(end); 
        trialID = [subjID '_' num2str(i)];
        
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
            sampleRayleighMatch([trialID '_threshold'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,testSpd,...
            'threshold','plotResults',...
            true,'nObserverMatches',p.Results.nObserverMatches,...
            'nReversals',p.Results.nReversals,'thresholdScaleFactor',...
            p.Results.thresholdScaleFactor,'nBelowThreshold',...
            p.Results.nBelowThreshold,'age',p.Results.age,'p1Scale',...
            p.Results.p1Scale,'p2Scale',p.Results.p2Scale,'testScale',...
            p.Results.testScale,'LMEqualOD',p.Results.LMEqualOD,...
            'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength);
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,testSpd,...
            'forcedChoice','plotResults',true,'nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength);
    end
    
    % Vary number of matches
elseif strcmp(testingParamToVary,'nMatches')
    paramName = 'Number of Matches Per Test Light';
    for i = 1:length(testingValsToVary)
        trialID = [subjID '_' num2str(i)];
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
            sampleRayleighMatch([trialID '_threshold'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,test,...
            'threshold','plotResults',true,'nObserverMatches',...
            testingValsToVary(i),'nReversals',p.Results.nReversals,...
            'thresholdScaleFactor',p.Results.thresholdScaleFactor,...
            'nBelowThreshold',p.Results.nBelowThreshold,'age',...
            p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength);
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,test,...
            'forcedChoice','plotResults',true,'nObserverMatches',...
            testingValsToVary(i),'nReversals',p.Results.nReversals,'age',...
            p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',p.Results.adjustmentLength);
    end
    % Vary simulation adjustment length. Note that this only runs the 
    % forced choice version, as threshold often cannot be reached for 
    % coarser arrays
elseif strcmp(testingParamToVary,'adjustmentLength')
    paramName = 'Number of Possible Light Adjustments';
    for i = 1:length(testingValsToVary)
        trialID = [subjID '_' num2str(i)];
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch([trialID '_FC'],nObservers,...
            p.Results.baseConeParams,coneParamsToVary,p1,p2,test,...
            'forcedChoice','plotResults',true,'nObserverMatches',...
            p.Results.nObserverMatches,'nReversals',p.Results.nReversals,...
            'age',p.Results.age,'p1Scale',p.Results.p1Scale,'p2Scale',...
            p.Results.p2Scale,'testScale',p.Results.testScale,'LMEqualOD',...
            p.Results.LMEqualOD,'dlens0',p.Results.dlens0,'restrictBySd',...
            p.Results.restrictBySd,'dmac0',p.Results.dmac0,...
            'adjustmentLength',testingValsToVary(i));
    end
else
    error('Unrecognized parameter to vary');
end

% Save data 
save(fullfile(outputDir,'testData.mat'),'p','coneErrFC','matchErrFC',...
    'coneErrAdjust','matchErrAdjust','subjID','nObservers','p1','p2',...
    'test','coneParamsToVary','testingParamToVary','testingValsToVary');

%% Plotting
theFig = figure(1);
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
xVals = testingValsToVary;
% If wavelength increment is being varied, the x value is the number of
% wavelengths, not the increment.
if strcmp(testingParamToVary,'testWlIncr')
    wlRange = length(test(1):test(end));
    xVals = ceil(wlRange*(1./testingValsToVary));
end
% Plot results 
plot(xVals,coneErrFC,'r-o','LineWidth',2.5);
if ~strcmp(testingParamToVary,'adjustmentLength')
    plot(xVals,coneErrAdjust,'b-o','LineWidth',2.5);
    legend('Forced Choice','Adjustment');
end 
title('Average Cone Spectral Sensitivity Error');
xlabel(paramName);
ylabel('Average RMS Error');
ylim([0 0.02]);
text(xVals(end)/10,0.0175,['n = ' num2str(nObservers)]);

% Subplot 2 - Match error
subplot('Position', subplotPosVectors(2,1).v);
hold on;
plot(xVals,matchErrFC,'r-o','LineWidth',2.5);
if ~strcmp(testingParamToVary,'adjustmentLength')
    plot(xVals,matchErrAdjust,'b-o','LineWidth',2.5);
    legend('Forced Choice','Adjustment');
end 
title('Average Match Error');
xlabel(paramName);
ylabel('Error');
ylim([0 0.05]);
text(xVals(end)/10,0.045,['n = ' num2str(nObservers)]);

% Title and save plots 
sgtitle(['Vary ' paramName]);
NicePlot.exportFigToPDF(fullfile(outputDir,'plots'),theFig,300);
end