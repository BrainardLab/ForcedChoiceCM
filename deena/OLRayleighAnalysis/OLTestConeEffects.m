function OLTestConeEffects(fName, varargin)
% Calculate cone effects of subjects' OL Rayleigh matches.
%
% Syntax:
%   OLTestConeEffects(fName, varargin)
%
% Description
%    Takes in a file of user's Rayleigh matches on the OneLight (produced 
%    by OLRayleighMatch). For each match, calculates and plots a comparison 
%    of nominal and predicted spds, a comparison of test and match 
%    predicted spds, and a bar graph of cone responses to the test and 
%    match lights. Can also run to plot test/match cone responses based on 
%    the radiometer playback of subject matches. Figures are saved in a 
%    subject and session-specific folder.
%
% Inputs:
%    fName      - character array of filename. Ends in .mat
%
% Outputs:
%    Saves three figures for each match in the match file. When a radiometer
%    playback file is used, saves one ffigure for each match.  
%
% Optional key-value pairs:
%    'fType'     - character array of file type. Default is '.tif'
%    'measured'  - logical indicating whether to calculate cone effects
%                  from radiometer measurements (true) or nominal spd data
%                  (false). Default is false.

% Example: OLTestConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/Deena/Deena_1.mat')
% Example: OLTestConeEffects('/Users1/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/Deena/Deena_1.mat')
% Example: OLTestConeEffects('/Users1/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/Deena/Deena_1.mat','measured',true)


%% Parse input
close all; 
p = inputParser;
p.addParameter('fType', 'tiff', @(x) (ischar(x)));
p.addParameter('measured', false, @(x) (islogical(x)));
p.parse(varargin{:});
fType = p.Results.fType;
measured = p.Results.measured;

%% Load data file and radiometer file if appropriate
theData = load(fName);
if (isfield(theData,'p1') == 0 || isfield(theData,'p2') == 0  ||...
        isfield(theData,'test') == 0 || isfield(theData,'matches') == 0 ...
        || isfield(theData,'matchPositions') == 0 || isfield(theData, 'cal') == 0 ||...
        isfield(theData, 'primarySpdsNominal') == 0 || isfield(theData,...
        'primarySpdsPredicted') == 0 || isfield(theData,...
        'testSpdsNominal') == 0|| isfield(theData, 'testSpdsPredicted') == 0 ...
        || isfield(theData,'primaryStartStops') == 0 || isfield(theData,...
        'testStartStops') == 0 || isfield(theData,'subjectID') == 0 || ...
        isfield(theData, 'sessionNum') == 0 || isfield(theData, 'whitePrimaries' == 0)...
        || isfield(theData, 'whiteSettings' == 0) || isfield(theData, 'whiteStarts' == 0)...
        || isfield(theData, 'whiteStops' == 0) || isfield(theData, 'whiteSpdNominal' == 0) ...
        || isfield(theData, 'annulusData' == 0))
    error('Data file does not contain required variables');
end

if measured
    fprintf('\n******** Loading radiometer file ********\n');
    measFile = [fName(1:end - 4), '_meas', '.mat'];
    if ~exist(measFile, 'file')
        error('Radiometer measurements not available'); 
    end 
    measData = load(measFile);
    if (isfield(measData,'measuredTestSpds') == 0 ||...
            isfield(measData,'measuredPrimarySpds') == 0 ||...
            isfield(measData,'measuredWhite') == 0)
        error('Radiometer file does not contain required variables');
    end
    fprintf('Radiometer measurements successfully loaded\n');
end

%% Make directory for saving files
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'), 'coneResponsePlots', theData.subjectID, num2str(theData.sessionNum));
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

%% Initialize arrays
% Calculate sizes and initialize spd arrays
[nMatches, ~] = size(theData.matches);
wls = theData.cal.computed.pr650Wls;
S = WlsToS(wls);
inc = wls(2) - wls(1);

% Generate standard cone fundamentals for observer
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc);

% Initialize arrays
primaryCones = zeros(3, nMatches);
testCones = zeros(3, nMatches);

%% Calculate and plot effects of spectra on cones
for i = 1:nMatches
    if ~measured
        % Calculate effects of the spectra on cones
        testSpdNominal = theData.testSpdsNominal(:,theData.matchPositions(i,1));
        testSpdPredicted = theData.testSpdsPredicted(:,theData.matchPositions(i,1));
        testCones(:,i) = T_cones * testSpdPredicted;
        
        primarySpdNominal = theData.primarySpdsNominal(:,theData.matchPositions(i,2));
        primarySpdPredicted = theData.primarySpdsPredicted(:,theData.matchPositions(i,2));
        primaryCones(:,i) = T_cones * primarySpdPredicted;
        
        % Plots comparing nominal and predicted spds
        figure;
        subplot(1,2,1); hold on
        plot(wls,testSpdPredicted,'r','LineWidth',4);
        plot(wls,testSpdNominal,'k','LineWidth',2);
        title('Test');
        
        subplot(1,2,2); hold on
        plot(wls,primarySpdPredicted,'r','LineWidth',4);
        plot(wls,primarySpdNominal,'k','LineWidth',2);
        title('Primaries');
        legend('Predicted', 'Nominal');
        theTitle = sprintf('Match %g Nominal and Predicted Spds', i);
        % sgtitle(theTitle);
        file = fullfile(outputDir, [theData.subjectID, '_', theTitle]);
        saveas(gcf, file, fType);
        
        % Plots comparing test and match spds
        figure();
        plot(wls, testSpdPredicted, 'r', wls, primarySpdPredicted, 'b');
        theTitle = sprintf('Match %g Nominal Spds', i);
        title(theTitle);
        legend({ 'test' 'primaries' });
        file = fullfile(outputDir, [theData.subjectID, '_', theTitle]);
        saveas(gcf, file, fType);
    else
        testCones(:,i) = T_cones * measData.measuredTestSpds(:, i);
        primaryCones(:,i) = T_cones * measData.measuredPrimarySpds(:, i);
    end
    
    % Plot cone effects
    figure();
    cones = [testCones(1,i), primaryCones(1,i); testCones(2,i),...
        primaryCones(2,i); testCones(3,i), primaryCones(3,i)];
    bar(cones);
    
    names ={'L'; 'M'; 'S' };
    set(gca,'xticklabel', names)
    if measured 
        category = 'Measured'; 
    else 
        category = 'Nominal'; 
    end 
    theTitle = sprintf('Match %g %s Cone Responses', i, category); 
    title(theTitle);
    legend('Test','Primaries');
    file = fullfile(outputDir, [theData.subjectID, '_', theTitle]);
    saveas(gcf, file, fType);
end
end