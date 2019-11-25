% Script for computing cone responses of OL matches. Used to verify whether
% subjects' settings for primary ratio and test intensity in a given match
% lead to similar effects on their cones.

% Ask user to enter filename, load data
fName = input('Enter match filename: ');
load(fName);
if (exist('p1', 'var') == 0 || exist('p2', 'var') == 0 ||...
        exist('test', 'var') == 0 || exist('matches', 'var') == 0)
    error('Passed file does not contain required variables');
end

% Calculate sizes and initialize spd arrays
[row, col] = size(matches);
wls = cal.computed.pr650Wls;
inc = wls(2) - wls(1);
fullWidthHalfMax = 20;

% Generate standard cone fundamentals for observer
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc);

% Initialize arrays
primariesSpd = zeros(length(wls), row);
testSpd = zeros(length(wls), row);

primariesCones = zeros(3, row);
testCones = zeros(3, row);

for i = 1:row
    % For each match, add primaries and test intensities to spectrum
    p1Spd = OLMakeMonochromaticSpd(cal, p1, fullWidthHalfMax);
    p2Spd = OLMakeMonochromaticSpd(cal, p2, fullWidthHalfMax);
    primariesSpd(:,i) = (matches(i, 2) * p1Spd) +...
        ((1 - matches(i, 2)) * p2Spd);
    
    testSpd(:,i) = OLMakeMonochromaticSpd(cal, test, fullWidthHalfMax)...
        * matches(i,1);
    
    % Calculate effects of the spectra on cones
    primariesCones(:,i) = T_cones * primariesSpd;
    testCones(:,i) = T_cones * testSpd;
    
    % Plot spds
    figure();
    OLplotSpdCheck(testSpd(:, i),cal);
    hold on;
    OLplotSpdCheck(primariesSpd(:, i), cal);
    theTitle = sprintf('Match %g Spds', i);
    title(theTitle);
    legend({ 'test' 'primaries' });
    
    % Plot cone effects
    figure();
    cones = [testCones(1,i), primariesCones(1,i); testCones(2,i), primariesCones(2,i); testCones(3,i), primariesCones(3,i)];
    bar(cones);
    
    somenames={'L'; 'M'; 'S' };
    set(gca,'xticklabel',somenames)
    theTitle = sprintf('Match %g Cone Responses', i);
    title(theTitle);
    legend('Test','Primaries');
end