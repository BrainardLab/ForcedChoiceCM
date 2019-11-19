% Ask user to enter filename, load data
fName = input('Enter match filename: ');
load(fName);
if (exist('p1', 'var') == 0 || exist('p2', 'var') == 0 ||...
        exist('test', 'var') == 0 || exist('matches', 'var') == 0)
    error('Passed file does not contain required variables'); 
end 

% Get calibration information 
cal = OLGetCalibrationStructure;


% Generate standard cone fundamentals for observer
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments);

% Calculate sizes and initialize spd arrays  
[row, col] = size(matches); 
S = [380 1 401];
wls = SToWls(S);
fullWidthHalfMax = 20;

testIndex = find(wls == test);
p1Index = find(wls == p1);
p2Index = find(wls == p2);

primariesSpectrum = zeros(length(wls), row); 
primariesCones = zeros(3, row); 
testSpectrum = zeros(length(wls), row);
testCones = zeros(3, row); 


for i = 1:row 
    % For each match, add primaries and test intensities to spectrum 
    primariesSpectrum(p1Index,i) = matches(i,2);
    primariesSpectrum(p2Index, i) = 1 - matches(i, 2); 
    testSpectrum(testIndex,i) = matches(i,1); 
    
    % Calculate effects of the spectra on cones
    primariesCones(:,i) = T_cones * primariesSpectrum; 
    testCones(:,i) = T_cones * testSpectrum; 
    
    % Plot results 
    figure(); 
    plot(wls, primariesCones(:,i), 'r', wls, testCones(:,i), 'b'); 
end 


    
    % load OLSampleMatches11_14_DB.mat;
    % DBMatches = matches;
    % load OLSampleMatches11_14_DE.mat;
    % DEMatches = matches;