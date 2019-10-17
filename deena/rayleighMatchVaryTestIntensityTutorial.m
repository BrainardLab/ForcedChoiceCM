% Visualization of how varying test wavelength can help distinguish between
% observers who make similar Rayleigh matches but have different underlying
% lambda maxes and optical densities. Uses the plotRayleighMatches and 
% ComputeRayleighConfusions to produce a Pitt diagram of result. 

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
dphotopigments = [  ...
    [0 0 0]' ...
    [-90 -90 0]'];

colors = [ 'r' 'g' ];
legend = { 'observer1' 'observer2' };
title = 'Optical Density Variation';
plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments, colors,...
    legend, title);
