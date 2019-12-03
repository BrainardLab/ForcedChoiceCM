% Calculates cone effects of subjects' OL Rayleigh matches
function OLTestConeEffects(fName)
% Syntax:
%   OLTestConeEffects(fName)
%
% Description
%    Takes in a file of user's Rayleigh matches on the OneLight. For each \
%    match, the program calculates and plots a comparison of nominal and
%    predicted spds, a comparison of test and match predicted spds, and a
%    bar graph of cone responses to the test and match lights. Figures are
%    saved as PDFs in a subject-specific folder in MELA_analysis. 
%
% Inputs:
%    fName  - character array of filename. Ends in .mat
%
% Outputs:
%    saves three PDFs for each match in the match file
%
% Optional key-value pairs:
%    none 

% Example: OLTestConeEffefcts('/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatches/test/test_1.mat')

theData = load(fName);
if (isfield(theData,'p1') == 0 || isfield(theData,'p2') == 0  ||...
        isfield(theData,'test') == 0 || isfield(theData,'matches') == 0 ...
        || isfield(theData,'matchPositions') == 0 || isfield(theData, 'cal') == 0 ||...
        isfield(theData, 'primarySpdsNominal') == 0 || isfield(theData,...
        'primarySpdsPredicted') == 0 || isfield(theData,...
        'testSpdsNominal') == 0|| isfield(theData, 'testSpdsPredicted') == 0 ...
        || isfield(theData,'primaryStartStops') == 0 || isfield(theData,...
        'testStartStops') == 0 || isfield(theData,'subjectID') == 0 || ...
        isfield(theData, 'sessionNum') == 0)
    error('Data file does not contain required variables');
end

% Make directory for saving files 
outputDir = fullfile(getpref('ForcedChoiceCM','rayleighAnalysisDir'), 'cone response plots', theData.subjectID);
if (~exist(outputDir,'dir'))
    mkdir(outputDir);
end

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

for i = 1:nMatches
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
    theTitle = sprintf('Match %g nominal and predicted spds', i);
    % sgtitle(theTitle);
    file = fullfile(outputDir, [theData.subjectID, '_', theTitle, '.pdf']); 
    saveas(gcf, file); 
    
    % Plots comparing test and match spds
    figure();
    plot(wls, testSpdPredicted, 'r', wls, primarySpdPredicted, 'b');
    theTitle = sprintf('Match %g Spds', i);
    title(theTitle);
    legend({ 'test' 'primaries' });
    file = fullfile(outputDir, [theData.subjectID, '_', theTitle, '.pdf']); 
    saveas(gcf, file); 
    
    % Plot cone effects
    figure();
    cones = [testCones(1,i), primaryCones(1,i); testCones(2,i),...
        primaryCones(2,i); testCones(3,i), primaryCones(3,i)];
    bar(cones);
    
    names ={'L'; 'M'; 'S' };
    set(gca,'xticklabel', names)
    theTitle = sprintf('Match %g Cone Responses', i);
    title(theTitle);
    legend('Test','Primaries');
    file = fullfile(outputDir, [theData.subjectID, '_', theTitle, '.pdf']); 
    saveas(gcf, file); 
end
end