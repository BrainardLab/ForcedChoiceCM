% Load light settings (right now I am using this for measured matches, but
% can edit it to be for nominal)
fName = '/Users/melanopsin/Dropbox (Aguirre-Brainard Lab)/MELA_data/Experiments/ForcedChoiceCM/OLRayleighMatch/ideal_meas.mat'; 
lightSettings = load(fName);

nonMeasured = load(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', 'OLRayleighMatchFineSpectralSettings_670_560_600.mat'));

% Generate standard cone fundamentals for observer. We will start with the
% foveal match case
foveal = true; 
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
inc = 2;
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc,...
    'foveal', foveal);

% Initialize arrays
[primaryRow, primaryCol] = size(lightSettings.primaryData);
[testRow, testCol] = size(lightSettings.testData);

primaryConeEffects = zeros(primaryCol, 3);
testConeEffects = zeros(testCol, 3);


for i = 1:primaryCol
    primaryConeEffects(i,:) = (T_cones * lightSettings.primaryData(:,i))';
end

for i = 1:testCol
    testConeEffects(i,:) = (T_cones * lightSettings.testData(:,i))';
end 

% Find least squares error
minErr = 1e8;
tIndex = 0;
pIndex = 0;

for i = 1:testCol
    for j = 1:primaryCol
        err  = (testConeEffects(i,1) - primaryConeEffects(j,1) )^2 + ...
            (testConeEffects(i,2) - primaryConeEffects(j,2))^2 +...
            (testConeEffects(i,3) - primaryConeEffects(j,3)) ^2;
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
fprintf('Primary Setting: %g \n',nonMeasured.p1Scales(pIndex+180));

% Plot its cone effects and spd
figure;
OLPlotConeEffects(primaryConeEffects(pIndex,:)',testConeEffects(tIndex,:)','Ideal',1);
figure;
OLplotSpdCheck(380:2:780,[lightSettings.testData(:,tIndex)...
    lightSettings.primaryData(:,pIndex)]);
legend('test','primary');
title('Measured Spds for Ideal Match');

% Display match on OneLight
ol = OneLight;
for i = 1:10
    ol.setMirrors(squeeze(nonMeasured.primaryStartStops(pIndex+180,1,:))',...
        squeeze(nonMeasured.primaryStartStops(pIndex+180,2,:))');
    pause(0.5);
    ol.setMirrors(squeeze(nonMeasured.testStartStops(tIndex,1,:))',...
        squeeze(nonMeasured.testStartStops(tIndex,2,:))');
    pause(0.5); 
end