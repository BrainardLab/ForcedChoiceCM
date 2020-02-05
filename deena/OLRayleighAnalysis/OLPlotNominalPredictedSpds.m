function OLPlotNominalPredictedSpds(fName)
% Plot nominal and predicted spds for subjects' OL Rayleigh matches.
%
% Syntax:
%   OLPlotNominalPredictedSpds(fName)
%
% Description
%    Plots nominal and predicted spds for subjects' OneLight Rayleigh 
%    matches, as recorded in their data files from the experiment. Produces
%    two plots for each match-one comparing nominal and predicted spds, and
%    one comparing predicted test and match spds. 
%
% Inputs:
%    fName      - Character array of filename produced by OLRayleighMatch. 
%                 Ends in .mat
%
% Outputs:
%    Two plots per match, as described above. 
%
% Optional key-value pairs:
%    none 

% History 
%    1/22/20   dce  -Modified program from OLTestConeEffects   
% Example: OLGetConeEffects('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatches/test/test_1.mat')

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

% Make directory for saving files 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'), 'coneResponsePlots', theData.subjectID, num2str(theData.sessionNum));
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

% Extract useful parameters from data
[nMatches, ~] = size(theData.matches);
wls = theData.cal.computed.pr650Wls;

% Make plots for each match
for i = 1:nMatches
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
    % sgtitle(theTitle);
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
end
