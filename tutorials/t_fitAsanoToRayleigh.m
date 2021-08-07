% t_FitAsanoToRayleigh
%
% Description:
%    Demonstration/test tutorial to check/develop fitting methods to
%    Rayleigh match data.
%

% History:
%    07/29/21  dhb  Started.

%% Clear
clear; close all;

%% Load and homogenize data
whichData = 'simulatedNoNoise';
theDir = fileparts(mfilename('fullpath'));
switch (whichData)
    case 'simulatedNoise'
        dataFilename = 'SimulatedWithNoise_spds.mat';
        theData = load(fullfile(theDir,'sampleData',dataFilename));
        S = [380 2 201]; wls = SToWls(S);
        fieldSize = theData.fieldSize;
        age = theData.age;
        [nWls,nReps,nRefs] = size(theData.primarySpdsByWl);
        refSpds = reshape(theData.refSpdsByWl,S(3),nReps*nRefs);
        primarySpds = reshape(theData.primarySpdsByWl,S(3),nReps*nRefs);
        standardParams = zeros(1,8);
        simulatedParams = theData.sampledParams;
    case 'simulatedNoNoise'
        dataFilename = 'SimulatedWithoutNoise_3201Levels.mat';
        %dataFilename = 'SimulatedWithoutNoise_noiselessSearch.mat';
        theData = load(fullfile(theDir,'sampleData',dataFilename));
        S = [380 2 201]; wls = SToWls(S);
        fieldSize = theData.fieldSize;
        age = theData.age;
        [nWls,nReps,nRefs] = size(theData.primarySpdsByWl);
        refSpds = reshape(theData.refSpdsByWl,S(3),nReps*nRefs);
        primarySpds = reshape(theData.primarySpdsByWl,S(3),nReps*nRefs);
        standardParams = zeros(1,8);
        simulatedParams = theData.sampledParams;
    case 'exactMatches'
        dataFilename = 'SimulatedWithoutNoise_spds.mat';
        theData = load(fullfile(theDir,'sampleData',dataFilename));
        S = [380 2 201]; wls = SToWls(S);
        fieldSize = theData.fieldSize;
        age = theData.age;
        [nWls,nReps,nRefs] = size(theData.primarySpdsByWl);
        refSpds = reshape(theData.refSpdsByWl,S(3),nReps*nRefs);
        primarySpds = reshape(theData.primarySpdsByWl,S(3),nReps*nRefs);
        standardParams = zeros(1,8);
        simulatedParams = theData.sampledParams;

        % Use the data we read in to generate some apparatus primaries


    case 'MELC_0004'
        dataFilename = 'MELC_0004_spds.mat';
        theData = load(fullfile(theDir,'sampleData',dataFilename));
        S = [380 2 201]; wls = SToWls(S);
        fieldSize = 10;
        age = 32;
        [nWls,nReps,nRefs] = size(theData.measPrimarySpdsByWl);
        refSpds = reshape(theData.measRefSpdsByWl,S(3),nReps*nRefs);
        primarySpds = reshape(theData.measPrimarySpdsByWl,S(3),nReps*nRefs);
        standardParams = zeros(1,8);
        simulatedParams = zeros(1,8);
    otherwise
        error('What is this data file of which you speak?')
end

%% Standard and simulated observer I
standardObserver = genRayleighObserver('fieldSize',fieldSize,'age',...
    age,'calcCones',true,'coneVec', standardParams,...
    'S',S,'opponentParams',[40.3908 205.7353 62.9590 1.0000]);
simulatedObserver = genRayleighObserver('fieldSize',fieldSize,'age',...
    age,'calcCones',true,'coneVec', simulatedParams,...
    'S',S,'opponentParams',[40.3908 205.7353 62.9590 1.0000]);

%% Call Deena fit routine with robust fitting loop
minimizeConeError = false;
maxIterations = 1;
trimSd = 1;
fitRefSpds{1} = refSpds;
fitPrimarySpds{1} = primarySpds;
for jj = 1:maxIterations
    nDataLeftToFit = size(fitRefSpds{jj},2);
    fitParamsIter{jj} = findObserverParameters(fitRefSpds{jj},fitPrimarySpds{jj},'S0',true,'dmac0',true,'sdLambdaMax',4, ...
        'age',age,'fieldSize',fieldSize, ...
        'minimizeConeErr',false,'matchErrorFun',@findMatchErrorNew,'initialConeParams',simulatedParams);
    [fitObjective{jj},fitObserverErrors{jj}] = findMatchErrorNew(fitParamsIter{jj},standardObserver,fitRefSpds{jj},fitPrimarySpds{jj}, ...
        'findConeErr',minimizeConeError);
    errorStdevs{jj} = std(fitObserverErrors{jj},0,2);
    includeIndex = find(fitObserverErrors{jj}(1,:) <= trimSd*errorStdevs{jj}(1) & fitObserverErrors{jj}(2,:) <= trimSd*errorStdevs{jj}(2));
    if (length(includeIndex) == nDataLeftToFit)
        break;
    end
    if (isempty(includeIndex))
        error('Oops, no data left to fit')
    end
    fitRefSpds{jj+1} = fitRefSpds{jj}(:,includeIndex);
    fitPrimarySpds{jj+1} = fitPrimarySpds{jj}(:,includeIndex);
end
fitParams = fitParamsIter{end};

% Get observer with the final parameters
fitObserver = ObserverVecToParams('basic', ...
    [fitParams standardObserver.colorDiffParams.noiseSd],standardObserver);
fitObserver.T_cones = ComputeObserverFundamentals(fitObserver.coneParams,...
    S);

%% Calculate and plot differences
[standardObjective,standardObserverErrors] = findMatchErrorNew(standardParams,standardObserver,refSpds,primarySpds, ...
    'findConeErr',minimizeConeError);
[simulatedObjective,simulatedObserverErrors] = findMatchErrorNew(simulatedParams,simulatedObserver,refSpds,primarySpds, ...
    'findConeErr',minimizeConeError);
[fitObjective,fitObserverErrors] = findMatchErrorNew(fitParams,standardObserver,refSpds,primarySpds, ...
    'findConeErr',minimizeConeError);

figure; clf;
if (minimizeConeError)
    axisLim = 0.05;
else
    axisLim = 15;
end
subplot(1,3,1); hold on
index = 1;
colorIndex = 1;
theColors = ['r' 'g' 'b' 'k' 'c' 'm' 'y' ];
for ii = 1:nRefs
    for jj = 1:nReps
        plot(standardObserverErrors(1,index),standardObserverErrors(2,index),[theColors(colorIndex) 'o'],'MarkerFaceColor',theColors(colorIndex),'MarkerSize',10);
        index = index+1;
        colorIndex = colorIndex + 1;
        if (colorIndex > length(theColors))
            colorIndex = 1;
        end
    end
end
axis('square');
plot([-axisLim axisLim],[0 0],'k:','LineWidth',0.5);
plot([0 0],[-axisLim axisLim],'k:','LineWidth',0.5);
xlim([-axisLim axisLim]); ylim([-axisLim axisLim]);
title('Standard')

subplot(1,3,2); hold on
index = 1;
colorIndex = 1;
theColors = ['r' 'g' 'b' 'k' 'c' 'm' 'y' ];
for ii = 1:nRefs
    for jj = 1:nReps
        plot(simulatedObserverErrors(1,index),simulatedObserverErrors(2,index),[theColors(colorIndex) 'o'],'MarkerFaceColor',theColors(colorIndex),'MarkerSize',10);
        index = index+1;
        colorIndex = colorIndex + 1;
        if (colorIndex > length(theColors))
            colorIndex = 1;
        end
    end
end
axis('square');
plot([-axisLim axisLim],[0 0],'k:','LineWidth',0.5);
plot([0 0],[-axisLim axisLim],'k:','LineWidth',0.5);
xlim([-axisLim axisLim]); ylim([-axisLim axisLim]);
title('Simulated')

subplot(1,3,3); hold on
index = 1;
colorIndex = 1;
theColors = ['r' 'g' 'b' 'k' 'c' 'm' 'y' ];
for ii = 1:nRefs
    for jj = 1:nReps
        plot(fitObserverErrors(1,index),fitObserverErrors(2,index),[theColors(colorIndex) 'o'],'MarkerFaceColor',theColors(colorIndex),'MarkerSize',10);
        index = index+1;
        colorIndex = colorIndex + 1;
        if (colorIndex > length(theColors))
            colorIndex = 1;
        end
    end
end
axis('square');
plot([-axisLim axisLim],[0 0],'k:','LineWidth',0.5);
plot([0 0],[-axisLim axisLim],'k:','LineWidth',0.5);
xlim([-axisLim axisLim]); ylim([-axisLim axisLim]);
title('Fit')

function [errorToMinimize,errors] = findMatchErrorNew(coneParamsVec,initialObs,testSpds,primarySpds,varargin)
% Computes the error associated with a set of Rayleigh matches
%
% Syntax:
%   [errorToMinimize,errors] = findMatchErrorNew(coneParamsVec,initialObs,testSpds,primarySpds);
%
% Description:
%    Given observer cone parameters, computes the error associated with
%    a set of pairs of spectra which the observer identified as Rayleigh
%    matches. For each pair, computes the observer's cone responses to both
%    test and primary lights, then converts this to opponent contrast. The
%    error for a given pair is represented as the vector length of the
%    luminance and RG opponent contrast terms. The overall error is the 
%    root mean square of the individual error terms.
%
%    This function is designed for use with fmincon or similar parameter
%    search functions.
%
% Inputs:
%    coneParamsVec -Vector of eight individual difference parameters. See
%                 ObserverVecToParams for a full description
%    initialObs  -Struct containing the initial settings for the observer.
%                 Not modified by the program, but some fields are used
%                 for reference.
%    testSpd     -Vector representation of the predicted spds for
%                 the chosen test light.
%    primarySpd  -Vector representation of the predicted spds for
%                 the chosen primary light.
%
% Outputs:
%    error       -Root mean square of the vector lengths of luminance and
%                 RG opponent contrast terms.
%
% Optional key-value pairs:
%    S           -Wavelength sampling for cone calculations, in the 
%                 form [start increment numTerms]. Default is [380 2 201] 
%                 (OneLight convention)  
%    errScalar   -Integer for scaling the match error, to improve search.
%                 Default is 100.
%    findConeErr -Logical. If true, calculates cone exictation error
%                 instead of opponent contrast difference. Default 
%                 is false.

% History:
%   06/12/20  dce       Wrote it.
%   06/15/20  dce       Modified to take in multiple spds
%   05/09/21  dce       Added option to find error based on cone excitation
%                       difference instead of opponent contrast.

% Parse input 
p = inputParser;
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('errScalar',100,@(x)(isnumeric(x)));
p.addParameter('findConeErr',false,@(x)(islogical(x)));
p.parse(varargin{:});

% Throw error if matrix sizes do not match 
[spdLength,nMatches] = size(testSpds);   
if length(SToWls(p.Results.S)) ~= spdLength
    error('Observer wavelength sampling and spd length do not match'); 
end 

% Array for storing error of each pair
errors = zeros(3,nMatches);
pairError = zeros(1,nMatches);   

% Find opponent parameters
CDParams = initialObs.colorDiffParams;
opponentParams = [CDParams.lumWeight CDParams.rgWeight CDParams.byWeight...
    CDParams.noiseSd];

% Generate an observer
observer = genRayleighObserver('age', initialObs.coneParams.ageYears,...
    'fieldSize', initialObs.coneParams.fieldSizeDegrees,...
    'coneVec',coneParamsVec,'opponentParams',opponentParams,'S',p.Results.S); 

% Calculate error metric for each match 
for i = 1:nMatches
    % Calculate cone responses for the given spectra
    test_LMS = observer.T_cones * testSpds(:,i);
    primary_LMS = observer.T_cones * primarySpds(:,i);
    
    if p.Results.findConeErr
        % Calculate cone contrast difference
        errors(:,i) = 0.5*(test_LMS - primary_LMS) ./ (test_LMS + primary_LMS);
    else
        % Calculate opponent contrast
        opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
            test_LMS, primary_LMS);
        errors(:,i) = opponentContrast;
    end

    % Find vector length of error (excluding S components)
    pairError(i) = norm(errors(1:2,i));
end

% Report the root mean square error (scaled up to improve searching)
errorToMinimize = sqrt(mean(pairError.^2))*p.Results.errScalar; 

end