% Script for finding the "ideal" match on the OneLight for a standard
% observer, given the measured primary lights used in the experiment.

% History
%    dce    xx/xx/20  - Wrote it
%    dce    4/5/20    - Edited for style
%    dce    6/9/20    - Modified to work with calculated spds
%    dce    6/16/20   - Modified to use simulated observer function

% Load light settings - currently set to measured data.
% fName = '/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/ideal_meas.mat';
% lightSettings = load(fName);

% Uncomment these lines to compute ideal match based on nominal rather than
% measured data (may need to change some other variable names).
nonMeasured = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops',...
    'OLRayleighMatchFineSpectralSettings_670_560_600_1_0.02_0.07.mat'));
lightSettings = nonMeasured;

% Generate cone fundamentals for observer. Can adjust parameters
% as needed.
age = 32;
fieldSizeDeg = 2;
coneParams = [0 0 0 0 0 4 0 0 0];
observer = genRayleighObserver('age',age,'fieldSize',fieldSizeDeg,'coneVec',coneParams);
T_cones = observer.T_cones;

% foveal = true;                          % Foveal matching
% lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
% dphotopigments = [0 0 0]';              % No optical density variation
% inc = 2;                                % Increment between wavelengths in spd measurements
wls = 380:2:780;                        % Range of wavelengths tested
% T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc,...
%     'foveal', foveal);

% Initialize arrays
[~, primaryCol] = size(lightSettings.primarySpdsPredicted);
[~, testCol] = size(lightSettings.testSpdsPredicted);

primaryConeEffects = zeros(primaryCol, 3);
testConeEffects = zeros(testCol, 3);

% Calculate cone effects for each primary/test light in the dataset
for i = 1:primaryCol
    primaryConeEffects(i,:) = (T_cones * lightSettings.primarySpdsPredicted(:,i))';
end

for i = 1:testCol
    testConeEffects(i,:) = (T_cones * lightSettings.testSpdsPredicted(:,i))';
end

%% Find least squares error

% Initial settings for minimum error and for primary/test indices
minErr = Inf;
tIndex = 0;
pIndex = 0;

% Loop through each primary/test light combination and calculate the
% squared error of the difference in cone effects
for i = 1:testCol
    for j = 1:primaryCol
        opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
            testConeEffects(i,:)', primaryConeEffects(j,:)');
        err = norm(opponentContrast(1:2)); 
        
        if err < minErr
            minErr = err;
            tIndex = i;
            pIndex = j;
        end
    end
end

% Show info aboutideal match
fprintf('Minimum Error: %g \n',minErr);
fprintf('Test Setting: %g \n',nonMeasured.testScales(tIndex));
fprintf('Primary Setting: %g \n',nonMeasured.p1Scales(pIndex));

% Plot cone effects and spds at the minimum error settings
OLPlotConeEffects(primaryConeEffects(pIndex,:)',...
    testConeEffects(tIndex,:)','Ideal',1);
OLPlotSpdCheck(wls,[lightSettings.testSpdsPredicted(:,tIndex)...
    lightSettings.primarySpdsPredicted(:,pIndex)]);
legend('test','primary');
title('Measured Spds for Ideal Match');

% Display match on OneLight - primary and test lights will appear in
% half-second intervals for 10 rotations
% ol = OneLight();
% for i = 1:10
%     ol.setMirrors(squeeze(nonMeasured.primaryStartStops(pIndex+180,1,:))',...
%         squeeze(nonMeasured.primaryStartStops(pIndex+180,2,:))');
%     pause(0.5);
%     ol.setMirrors(squeeze(nonMeasured.testStartStops(tIndex,1,:))',...
%         squeeze(nonMeasured.testStartStops(tIndex,2,:))');
%     pause(0.5);
% end