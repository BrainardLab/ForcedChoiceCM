% Load light settings
fName = fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedStartStops', 'OLRayleighMatchFineSpectralSettings.mat');
lightSettings = load(fName);

% Generate standard cone fundamentals for observer. We will start with the
% foveal match case
lambdaMaxes = [558.9 530.3 420.7]';     % Normal trichromat
dphotopigments = [0 0 0]';
inc = 2;
T_cones = findConeFundamentals(lambdaMaxes, dphotopigments, 'inc', inc,...
    'foveal', true);

% Initialize arrays
[spdRow, spdCol] = size(lightSettings.testSpdsPredicted);
testConeEffects = zeros(spdCol, 3);
primaryConeEffects = zeros(spdCol, 3);


for i = 1:spdCol
    testConeEffects(i,:) = (T_cones * lightSettings.testSpdsPredicted(:,i))';
    primaryConeEffects(i,:) = (T_cones * lightSettings.primarySpdsPredicted(:,i))';
end

% Find least squares error
minErr = 1e8;
tIndex = 0;
pIndex = 0;

for i = 1:spdCol
    for j = 1:spdCol
        err  = (testConeEffects(i,1) - primaryConeEffects(j,1) )^2 + ...
            (testConeEffects(i,2) - primaryConeEffects(j,2))^2 + ...
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
fprintf('Test Setting: %g \n',lightSettings.testScales(tIndex));
fprintf('Primary Setting: %g \n',lightSettings.p1Scales(pIndex));

% Plot its cone effects and spd
figure;
OLPlotConeEffects(primaryConeEffects(pIndex,:)',testConeEffects(tIndex,:)','Ideal',1);
figure;
OLplotSpdCheck(380:2:780,[lightSettings.testSpdsPredicted(:,tIndex)...
    lightSettings.primarySpdsPredicted(:,pIndex)]);
legend('test','primary');
title('Predicted Spds for Ideal Match');

% Display match on OneLight
ol = OneLight;
for i = 1:10
    ol.setMirrors(squeeze(lightSettings.primaryStartStops(pIndex,1,:))',...
        squeeze(lightSettings.primaryStartStops(pIndex,2,:))');
    pause(1);
    ol.setMirrors(squeeze(lightSettings.testStartStops(pIndex,1,:))',...
        squeeze(lightSettings.testStartStops(pIndex,2,:))');
    pause(1); 
end