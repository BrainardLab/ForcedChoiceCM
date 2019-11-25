% Script for computing cone responses of OL matches. Used to verify whether
% subjects' settings for primary ratio and test intensity in a given match 
% lead to similar effects on their cones. 

% Ask user to enter filename, load data
fName = input('Enter match filename: ');
theData = load(fName);
if (isfield(theData,'p1') == 0 || isfield(theData,'p2') == 0  ||...
       isfield(theData,'test') == 0 || isfield(theData,'matches') == 0)
    error('Data file does not contain required variables');
end

% Get calibration information
cal = OLGetCalibrationStructure;

% Calculate sizes and initialize spd arrays
[nMatches, ~] = size(theData.matches);
wls = cal.computed.pr650Wls;
S = WlsToS(wls);
inc = wls(2) - wls(1);
fullWidthHalfMax = 20;
lambda = 0.001;

% Generate standard cone fundamentals for observer
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc);

% Initialize arrays
primariesSpdNominal = zeros(length(wls), nMatches);
testSpdNominal = zeros(length(wls), nMatches);
primariesSpdPredicted = zeros(length(wls), nMatches);
testSpdPredicted = zeros(length(wls), nMatches);

primariesCones = zeros(3, nMatches);
testCones = zeros(3, nMatches);

for i = 1:nMatches
    % For each match, add primaries and test intensities to spectrum
    p1Spd = OLMakeMonochromaticSpd(cal, theData.p1, fullWidthHalfMax);
    p2Spd = OLMakeMonochromaticSpd(cal, theData.p2, fullWidthHalfMax);
    primariesSpdNominal(:,i) = (matches(i, 2) * p1Spd) +...
        ((1 - matches(i, 2)) * p2Spd);
    [~,primariesTemp,primariesSpdPredicted(:,i)] = OLSpdToSettings(cal, primariesSpdNominal(:,i), 'lambda', lambda);

    
    testSpdNominal(:,i) = OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)...
        * matches(i,1);
    [~,~,testSpdPredicted(:,i)] = OLSpdToSettings(cal, testSpdNominal(:,i), 'lambda', lambda);

    % Calculate effects of the spectra on cones
    primariesCones(:,i) = T_cones * primariesSpdPredicted;
    testCones(:,i) = T_cones * testSpdPredicted;
    
    figure;
    subplot(1,2,1); hold on
    plot(wls,testSpdPredicted(:,i),'r','LineWidth',4);
    plot(wls,testSpdNominal(:,i),'k','LineWidth',2);
    subplot(1,2,2); hold on
    plot(wls,primariesSpdPredicted(:,i),'r','LineWidth',4);
    plot(wls,primariesSpdNominal(:,i),'k','LineWidth',2);
    
    % Plot spds
    figure();
    OLplotSpdCheck(testSpdNominal(:, i),cal);

    hold on;
    OLplotSpdCheck(primariesSpdNominal(:, i),cal);
    OLplotSpdCheck(testSpdPredicted(:, i),cal); 
    OLplotSpdCheck(primariesSpdPredicted(:, i),cal);
    theTitle = sprintf('Match %g Spds', i);
    title(theTitle);
    legend({ 'test' 'primaries' });
    
    % Plot cone effects
    figure();
    h = zeros(6); 
    h(1:3) = plot(wls, testCones(1,i) * T_cones(1,:), 'r', wls, testCones(2,i)...
         * T_cones(2,:), 'r', wls, testCones(3,i) * T_cones(3,:), 'r');
    hold on; 
    h(4:6) = plot(wls, primariesCones(1,i) * T_cones(1,:), 'b', wls, primariesCones(2,i)...
        * T_cones(2,:), 'b', wls, primariesCones(3,i) * T_cones(3,:), 'b');
    theTitle = sprintf('Match %g Cone Responses', i);
    title(theTitle);
    legend(h(3:4),'Test','Primaries');
end