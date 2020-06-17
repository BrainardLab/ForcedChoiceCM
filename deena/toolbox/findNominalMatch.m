function [pSpd,tSpd,pIndex,tIndex] = findNominalMatch(fName,coneParams)
% Finds the nominal match for a Rayleigh matching experiment

% Syntax:
%   findNominalMatch(fName,coneParams)
%
% Description:
%    Takes in the file of light settings used to run the Rayleigh match 
%    experiment, as well as a nine-element vector of observer individual
%    difference parameters. Uses these inputs to compute which pair of 
%    primary and test lights has the smallest opponent contrast vector
%    length (which indicates that this pair is the best match). Returns
%    the primary and test spds, as well as their indices in the lights
%    array.
% Inputs:
%    fName       -Name of file containing Rayleigh match light settings
%    coneParams  -9-element vector containing Asano individual difference
%                 parameters
% Outputs:
%    pSpd        -Predicted primary spd for the match
%    tSpd        -Predicted test spd for the match
%    pIndex      -Index of the selected primary spd in the light settings
%                 array
%    tSpd        -Index of the selected test spd in the light settings
%                 array
%
% Optional key-value pairs:
%    None

% History
%    dce    xx/xx/20  - Wrote it
%    dce    4/5/20    - Edited for style
%    dce    6/9/20    - Modified to work with calculated spds
%    dce    6/16/20   - Modified to use simulated observer function

% Load light settings
lightSettings = load(fName);

% Generate cone fundamentals for observer. 
age = 32;
fieldSizeDeg = 2;
observer = genRayleighObserver('age',age,'fieldSize',...
    fieldSizeDeg,'coneVec',coneParams);
T_cones = observer.T_cones;

% Initialize arrays
[~, primaryCol] = size(lightSettings.primarySpdsPredicted);
[~, testCol] = size(lightSettings.testSpdsPredicted);

primaryConeEffects = zeros(primaryCol, 3);
testConeEffects = zeros(testCol, 3);

% Calculate cone effects for each primary/test light in the dataset
for i = 1:primaryCol
    primaryConeEffects(i,:) = (T_cones*lightSettings.primarySpdsPredicted(:,i))';
end

for i = 1:testCol
    testConeEffects(i,:) = (T_cones*lightSettings.testSpdsPredicted(:,i))';
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
fprintf('Test Setting: %g \n',lightSettings.testScales(tIndex));
fprintf('Primary Setting: %g \n',lightSettings.p1Scales(pIndex));

% Get spds of ideal match lights 
tSpd = lightSettings.testSpdsPredicted(:,tIndex); 
pSpd = lightSettings.testSpdsPredicted(:,pIndex); 

% Plot cone effects and spds at the minimum error settings
makePlots = true;
if makePlots
    OLPlotConeEffects(primaryConeEffects(pIndex,:)',...
        testConeEffects(tIndex,:)','Ideal',1);
    OLPlotSpdCheck(380:2:780,lightSettings.primarySpdsPredicted(:,pIndex));
    legend('test','primary');
    title('Measured Spds for Ideal Match');
end

end