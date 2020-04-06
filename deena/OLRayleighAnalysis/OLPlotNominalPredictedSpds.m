function OLPlotNominalPredictedSpds(fName, varargin)
% Plot nominal and predicted spds for subjects' OL Rayleigh matches.
%
% Syntax:
%   OLPlotNominalPredictedSpds(fName, varargin)
%
% Description
%    Plots nominal and predicted spds for subjects' OneLight Rayleigh
%    matches, as recorded in their data files from the experiment. Produces
%    two plots for each match-one comparing nominal and predicted spds, and
%    one comparing predicted test and match spds. When the option
%    'measured' is selected, makes two additional plots - one comparing
%    predicted and measured spds, and one comparing measured test and match
%    spds.
%
% Inputs:
%    fName       - Character array of filename produced by OLRayleighMatch.
%                  Ends in .mat
%
% Outputs:
%    Two plots per match, as described above.
%
% Optional key-value pairs:
%    'measured'  - Logical indicating whether to make plots of measured
%                  data in addition to the ones with nominal data. Default
%                  is false.
%    'measOnly'  - Logical indicating to only make plots of measured data,
%                  skipping over nominal vs predicted nominal test vs
%                  primary. Default is false.

% History
%    1/22/20   dce  - Modified program from OLTestConeEffects
%    4/5/20    dce  - Edited for style, added options for measured data
%
% Example: OLPlotNominalPredictedSpds('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatches/test/test_1.mat')

% Parse input
p = inputParser;
p.addParameter('measured', false, @(x) (islogical(x)));
p.addParameter('measOnly', false, @(x) (islogical(x)));
p.parse(varargin{:});

% Load data and check for needed variables
theData = load(fName);
if (isfield(theData,'matches') == 0 ...
        || isfield(theData,'matchPositions') == 0 ||...
        isfield(theData, 'cal') == 0 ||...
        isfield(theData, 'primarySpdsNominal') == 0 ||...
        isfield(theData, 'primarySpdsPredicted') == 0 ||...
        isfield(theData, 'testSpdsNominal') == 0 ||...
        isfield(theData, 'testSpdsPredicted') == 0 ||...
        isfield(theData,'subjectID') == 0 ||...
        isfield(theData, 'sessionNum') == 0 )
    error('Data file does not contain required variables');
end

% Load measured data if it is in use
if p.Results.measured
    fprintf('\n******** Loading radiometer file ********\n');
    measFile = [fName(1:end - 4), '_meas', '.mat'];
    if ~exist(measFile, 'file')
        error('Radiometer measurements not available');
    end
    measData = load(measFile);
    if (isfield(measData,'measuredTestSpds') == 0 ||...
            isfield(measData,'measuredPrimarySpds') == 0)
        error('Radiometer file does not contain required variables');
    end
    fprintf('Radiometer measurements successfully loaded\n');
end

% Make directory for saving files, based on subject ID and session number
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'),...
    'coneResponsePlots', theData.subjectID, num2str(theData.sessionNum));
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Extract useful parameters from data
[nMatches, ~] = size(theData.matches);  % Number of matches
wls = theData.cal.computed.pr650Wls;    % Wavelengths tested

% Make plots for each match
for i = 1:nMatches
    
    % Nominal plots. Skip if the 'measOnly' option is selected and measured
    % data is available. 
    if ~(p.Results.measured && p.Results.measOnly)
        % Extract spds for matches
        testSpdNominal = theData.testSpdsNominal(:,theData.matchPositions(i,1));
        testSpdPredicted = theData.testSpdsPredicted(:,theData.matchPositions(i,1));
        primarySpdNominal = theData.primarySpdsNominal(:,theData.matchPositions(i,2));
        primarySpdPredicted = theData.primarySpdsPredicted(:,theData.matchPositions(i,2));
        
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
        sgtitle(theTitle);
        file = fullfile(outputDir, [theData.subjectID, '_', strrep(theTitle,' ', '_')]);
        saveas(gcf, file, 'tiff');
        
        % Plots comparing test and match spds
        figure();
        plot(wls, testSpdPredicted, 'r', wls, primarySpdPredicted, 'b');
        theTitle = sprintf('Match %g Nominal Spds', i);
        title(theTitle);
        legend({ 'test' 'primaries' });
        file = fullfile(outputDir, [theData.subjectID, '_', strrep(theTitle,' ', '_')]);
        saveas(gcf, file, 'tiff');
    end
    
    % Plots for measured data
    if p.Results.measured
        testSpdMeasured = measData.measuredTestSpds(:, i);
        primarySpdMeasured = measData.measuredPrimarySpds(:, i);
        
        % Measured vs predicted spds
        figure;
        subplot(1,2,1); hold on
        plot(wls,testSpdPredicted,'r','LineWidth',4);
        plot(wls,testSpdMeasuredl,'k','LineWidth',2);
        title('Test');
        
        subplot(1,2,2); hold on
        plot(wls,primarySpdPredicted,'r','LineWidth',4);
        plot(wls,primarySpdMeasured,'k','LineWidth',2);
        title('Primaries');
        
        legend('Predicted', 'Measured');
        theTitle = sprintf('Match %g Measured and Predicted Spds', i);
        file = fullfile(outputDir, [theData.subjectID, '_', strrep(theTitle,' ', '_')]);
        saveas(gcf, file, 'tiff');
        
        % Measured test vs match spds
        figure();
        plot(wls, testSpdMeasured, 'r', wls, primarySpdMeasured, 'b');
        theTitle = sprintf('Match %g Measured Spds', i);
        title(theTitle);
        legend({ 'test' 'primaries' });
        file = fullfile(outputDir, [theData.subjectID, '_', strrep(theTitle,' ', '_')]);
        saveas(gcf, file, 'tiff');
    end
end
end
