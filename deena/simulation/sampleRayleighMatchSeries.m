function sampleRayleighMatchSeries(subjID,nObservers,baseConeParams,...
    coneParamsToVary,p1,p2,test,method,testingParamToVary,valsToVary,varargin)

%% Initial Setup
% Input parsing
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('LMEqualOD',false,@(x)(islogical(x)));
p.addParameter('dlens0',false,@(x)(islogical(x)));
p.addParameter('restrictBySd',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Make error-storing arrays
coneErrFC = zeros(length(valsToVary),1);
coneErrAdjust = zeros(length(valsToVary),1);
matchErrFC = zeros(length(valsToVary),1);
matchErrAdjust = zeros(length(valsToVary),1);

%% Define the variable parameter, then run the simulation in a loop
% Vary noise
if (length(testingParamToVary)==length('noise'))...
        && all(testingParamToVary=='noise')
    for i = 1:length(valsToVary)
        % Make base params with the specified noise level
        trialConeParams = [baseConeParams(1:8),valsToVary(i)];
        [coneErrAdjust(i),matchErrAdjust(i)] = ...
            sampleRayleighMatch(subjID,nObservers,trialConeParams,...
            coneParamsToVary,p1,p2,test,'threshold','plotResults',false);
        [coneErrFC(i),matchErrFC(i)] = ...
            sampleRayleighMatch(subjID,nObservers,trialConeParams,...
            coneParamsToVary,p1,p2,test,'forcedChoice','plotResults',false);
    end 
    
    % Vary wavelength increment
elseif (length(testingParamToVary)==length('testWlIncr'))...
        && all(testingParamToVary=='testWlIncr')
    for i = 1:length(valsToVary)
        testSpd = test(1):valsToVary(i):test(end);
        [coneErr(i),matchErr(i)] = ...
            sampleRayleighMatch(subjID,nObservers,baseConeParams,...
            coneParamsToVary,p1,p2,testSpd,method,'plotResults',false);
    end
    
    % Vary number of matches 
elseif (length(testingParamToVary)==length('testNMatches'))...
        && all(testingParamToVary=='testNMatches')
else
    error('Unrecognized parameter to vary');
end

% Plotting
xVals = valsToVary;
% If wavelength increment is being varied, the x value is the number of
% wavelengths, not the increment.
if (length(testingParamToVary)==length('testWlIncr'))...
        && all(testingParamToVary=='testWlIncr')
    wlRange = length(test(1):test(end));
    xVals = wlRange*(1./valsToVary);
end

% Cone spectral sensitivity error
figure();
plot(xVals,coneErrFC,'r-o',xVals,coneErrAdjust,'b-o','LineThickness',2);
theTitle = ['Average Cone Spectral Sensitivity Error - Vary ',...
    testingParamToVary];
title(theTitle);
xlabel('Parameter Value');
ylabel('Error');
legend('Forced Choice', 'Adjustment');

% Match error
figure();
plot(xVals,matchErrFC,'r-o',xVals,matchErrAdjust,'b-o','LineThickness',2);
theTitle = ['Average Match RMS Error - Vary ', testingParamToVary];
title(theTitle);
xlabel('Parameter Value');
ylabel('Error');
legend('Forced Choice', 'Adjustment');

p.addParameter('nObserverMatches',1,@(x)(isnumeric(x)));
p.addParameter('nReversals',[1 4],@(x)(isnumeric(x)));
p.addParameter('nBelowThreshold',1,@(x)(isnumeric(x)));
p.addParameter('thresholdScaleFactor',0.5,@(x)(isnumeric(x)));
end