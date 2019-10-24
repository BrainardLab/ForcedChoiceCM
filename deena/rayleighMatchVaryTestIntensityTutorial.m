% Visualization of how varying test wavelength can help distinguish between
% observers who make similar Rayleigh matches but have different underlying
% lambda maxes and optical densities. Uses the plotRayleighMatches and 
% ComputeRayleighConfusions to produce a Pitt diagram of result. 

%% Variation in cone optical density obscures differences in lambda max
%% for a hypothetical dichromat

% Cone lambda maxes, taken from Fig. 1 of the Thomas and Mollon paper.
lambdaMaxes = [  ...
    [541 541 420.7]' ...
    [551 551 420.7]'];

% We specify the cones as a shift relative to a nomogram generated lambda
%  max. These base values are set equal to the test values
baseLambdaMaxes = lambdaMaxes;

% You can also allow the specified photopigment density to vary. Enter
% these as percent changes relative to nominal values (positive or
% negative). Here, all photopigment densities are set to zero.
dphotopigments1 = [  ...
    [0 0 0]' ...
    [0 0 0]'];

dphotopigments2 = [  ...
    [0 0 0]' ...
    [-95 -95 0]'];
colors = [ 'r' 'g' ];
legend = { 'observer1: 541nm' 'observer2: 551nm' };

title = 'Optical Density Variation';
plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments1, colors,...
    legend, 'No Optical Density Variation');

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2, colors,...
    legend, 'Optical Density Variation: observer2 - 95%');

%% Can you separate the observers using different test lights?
plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2, colors,...
    legend, 'Optical Density Variation, test light to 570', 'testWavelength', 570);

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2, colors,...
    legend, 'Optical Density Variation, test light to 580', 'testWavelength', 580);

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2, colors,...
    legend, 'Optical Density Variation, test light to 600', 'testWavelength', 600);

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2, colors,...
    legend, 'Optical Density Variation, test light to 610', 'testWavelength', 610);

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2, colors,...
    legend, 'Optical Density Variation, test light to 620', 'testWavelength', 620);
