% Script for computing cone responses of OL matches. Used to verify whether
% subjects' settings for primary ratio and test intensity in a given match
% lead to similar effects on their cones.

% Ask user to enter filename, load data
fName = input('Enter match filename: ');
theData = load(fName);
if (isfield(theData,'p1') == 0 || isfield(theData,'p2') == 0  ||...
        isfield(theData,'test') == 0 || isfield(theData,'matches') == 0 ...
        || isfield(theData,'matchPositions') == 0 || isfield(theData, 'cal') == 0 ||...
        isfield(theData, 'primarySpdsNominal' == 0) || isfield(theData,...
        'primarySpdsPredicted' == 0) || isfield(theData,...
        'testSpdsNominal' == 0)|| isfield(theData, 'testSpdsPredicted' == 0) ...
        || isfield(theData,'primaryStartStops') == 0 || isfield(theData,'testStartStops') == 0)
    error('Data file does not contain required variables');
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
    theTitle = sprintf('Nominal and predicted spds: match %g', i);
    sgtitle(theTitle); 
    
    % Plots comparing test and match spds 
    figure();
    plot(wls, testSpdPredicted, 'r', wls, primarySpdPredicted, 'b'); 
    theTitle = sprintf('Match %g Spds', i);
    title(theTitle);
    legend({ 'test' 'primaries' });
    
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
end