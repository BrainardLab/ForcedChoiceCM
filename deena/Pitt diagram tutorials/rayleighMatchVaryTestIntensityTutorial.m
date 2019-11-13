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
theLegend = { 'observer1: 541nm' 'observer2: 551nm' };

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments1,...
    colors, theLegend, 'No Optical Density Variation, test light to 590nm');

plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2,...
    colors, theLegend, 'Optical Density Variation, test light to 590nm');

%% Plot the cone fundamentals for the two overlapping observers to see 
%  potential sites of modification

% Wavelength sampling. 
S = [380 1 401];
wls = SToWls(S); 

% Find cone fundamentals 
T_cones_observer1 = findConeFundamentals(lambdaMaxes(:,1),...
    dphotopigments1(:,1));
T_cones_observer2_noOD = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments1(:,2));
T_cones_observer2_OD = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments2(:,2));

figure; clf; hold on;
plot(wls,T_cones_observer1(1,:),'r','LineWidth',2);
plot(wls,T_cones_observer2_noOD(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (theLegend);
title ('Compare L cones - no OD variation'); 

figure; clf; hold on;
plot(wls,T_cones_observer1(1,:),'r','LineWidth',2);
plot(wls,T_cones_observer2_OD(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (theLegend);
title ('Compare L cones - OD variation'); 

% Find wavelength of maximum difference in each condition
[~, ind1] = max(abs(T_cones_observer1(1,:) - T_cones_observer2_noOD(2,:))); 
maxDiffWl1 = wls(ind1); 

[~, ind2] = max(abs(T_cones_observer1(1,:) - T_cones_observer2_OD(2,:))); 
maxDiffWl2 = wls(ind2); 


%% Use the wavelengths of max difference as test wavelengths to separate 
% the two observers 
plotTitle1 = sprintf('No Optical Density Variation, test light to %g nm'...
    , maxDiffWl1); 
plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments1,...
    colors, theLegend, plotTitle1, 'test', maxDiffWl1);

plotTitle2 = sprintf('Optical Density Variation, test light to %g nm',...
    maxDiffWl2); 
plotRayleighMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments2,...
    colors, theLegend, plotTitle2, 'test', maxDiffWl2);
