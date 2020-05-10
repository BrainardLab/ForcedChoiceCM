% Script for investigating how changing different lambda max/optical
% density parameters for trichromats affects simulated Rayleigh matches.
% Plots the ranges of these parameters across lambda max/optical density
% variation, then implements two case studies where optical density
% variations obscure lambda max variations. Identifies the wavelength with
% the largest difference in
% cone activations and uses this as a test light to separate the two cones.

%% Physiological ranges of parameters
% Cone lambda maxes - +-6nm for each cone
% Modified L cones
lambdaMaxesL = [  ...
    [552.9 530.3 420.7]'...
    [554.9 530.3 420.7]'...
    [556.9 530.3 420.7]'...
    [558.9 530.3 420.7]'... % Baseline
    [560.9 530.3 420.7]'...
    [562.9 530.3 420.7]'...
    [564.9 530.3 420.7]'];

lambdaMaxesL_small = [  ...
    [554.9 530.3 420.7]'...
    [558.9 530.3 420.7]'... % Baseline
    [562.9 530.3 420.7]'];

% Modified M cone
lambdaMaxesM = [  ...
    [558.9 524.3 420.7]'...
    [558.9 526.3 420.7]'...
    [558.9 528.3 420.7]'...
    [558.9 530.3 420.7]'... % Baseline
    [558.9 532.3 420.7]'...
    [558.9 534.3 420.7]'...
    [558.9 536.3 420.7]'];

lambdaMaxesM_small = [  ...
    [558.9 526.3 420.7]'...
    [558.9 530.3 420.7]'... % Baseline
    [558.9 534.3 420.7]'];

% Optical density variation (constrain to be the same for L and M)
ODs = [...
    [-30 -30 0]'...
    [-20 -20 0]' ...
    [-10 -10 0]'...
    [0 0 0]' ...
    [10 10 0]'...
    [20 20 0]' ...
    [30 30 0]'];

ODs_small = [...
    [-30 -30 0]'...
    [-15 -15 0]'...
    [0 0 0]' ...
    [15 15 0]'...
    [30 30 0]'];

[~, numODs] = size(ODs);
[~, numL] = size(lambdaMaxesL);
[~, numM] = size(lambdaMaxesM);
colors = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';

%% Case 1: Vary only lambda max (L and M separately)
% L cone variation
dphotopigments = repmat(ODs(:,4), 1, numL);
theTitle = sprintf('L cone variation, OD = 0 percent');
theLegend = {'552.9', '554.9', '556.9', '558.9', '560.9', '562.9', '564.9'};
plotRayleighMatches(lambdaMaxesL, lambdaMaxesL, dphotopigments,...
    colors(1:numL), theLegend, theTitle);

% M cone variation
dphotopigments = repmat(ODs(:,4), 1, numM);
theTitle = sprintf('M cone variation, OD = 0 percent');
theLegend ={'524.3', '526.3', '528.3', '530.3', '532.3', '534.3', '536.3'};
plotRayleighMatches(lambdaMaxesM, lambdaMaxesM, dphotopigments,...
    colors(1:numM), theLegend, theTitle);

%% Case 2: Vary both cones simultaneously (average OD)
% Include each cone's average lambda max and +-4nm
lambdaMaxes = [  ...
    [554.9 526.3 420.7]'...
    [558.9 526.3 420.7]'...
    [562.9 526.3 420.7]'...
    [554.9 530.3 420.7]'...
    [558.9 530.3 420.7]'... % Baseline
    [562.9 530.3 420.7]'...
    [554.9 534.3 420.7]'...
    [558.9 534.3 420.7]'...
    [562.9 534.3 420.7]'];

[~, n] = size(lambdaMaxes);
dphotopigments = repmat(ODs(:,4), 1, n); % Typical cones
theTitle = sprintf('Simultaneous cone variation, OD = 0 percent');
theLegend = {'554.9 +526.3', '558.9 +526.3', '562.9 +526.3', '554.9 +530.3',...
    '558.9 +530.3', '562.9 +530.3', '554.9 +534.3', '558.9 +534.3',...
    '562.9 +534.3'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments,...
    colors(1:n), theLegend, theTitle);

%% Case 3: Vary only OD (typical cones)
lambdaMaxes = repmat(lambdaMaxesL(:,4), 1, numODs); % typical cones
theTitle = sprintf('OD variation for standard cones');
theLegend = {'-30', '-20', '-10', '0', '10', '20', '30'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, ODs,...
    colors(1:numODs), theLegend, theTitle);

%% Case 4: Vary OD and one cone simultaneously
% M cone
lambdaMaxes = [];
dphotopigments = [];
for i = 1:length(lambdaMaxesM_small)
    lambdaMaxes = [lambdaMaxes, repmat(lambdaMaxesM_small(:,i), 1, length(ODs_small))];
    dphotopigments = [dphotopigments, ODs_small];
end
theLegend = {'526.3/-30', '526.3/-15', '526.3/0', '526.3/15', '526.3/30'...
    '530.3/-30', '530.3/-15', '530.3/0', '530.3/15', '530.3/30'...
    '534.3/-30', '534.3/-15', '534.3/0', '534.3/15', '534.3/30'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments,...
    colors, theLegend, 'Simultaneous Variation - M cone');

% L cone
lambdaMaxes = [];
dphotopigments = [];
for i = 1:length(lambdaMaxesL_small)
    lambdaMaxes = [lambdaMaxes, repmat(lambdaMaxesL_small(:,i), 1, length(ODs_small))];
    dphotopigments = [dphotopigments, ODs_small];
end
theLegend = {'554.9/-30', '554.9/-15', '554.9/0', '554.9/15', '554.9/30'...
    '558.9/-30', '558.9/-15', '558.9/0', '558.9/15', '558.9/30'...
    '562.9/-30', '562.9/-15', '562.9/0', '562.9/15', '562.9/30'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments,...
    colors, theLegend, 'Simultaneous Variation - L cone');

%% Examples
% Show that two distinguishable lambdaMaxes can be made to overlap or even
% reverse through OD variation, and try to separate them with test lights

% Wavelength sampling
S = [380 1 401];
wls = SToWls(S);
long_ind = 161:length(wls); % Indices of wavelengths above 540nm

%% Example 1 - M cone and OD variation
% lambdaMaxes +-2nm from average
lambdaMaxes = [
    [558.9 528.3 420.7]'...
    [558.9 532.3 420.7]'];
% Optical density +- 10%, 20% from average
dphotopigments1 = [...
    [0 0 0]'...
    [0 0 0]'];
dphotopigments2 = [...
    [10 10 0]'...
    [-10 -10 0]'];
dphotopigments3 = [...
    [20 20 0]'...
    [-20 -20 0]'];

% Legends
l1 = {'528.3', '532.3'};
l2 = {'528.3/+10', '532.3/-10'};
l3 = {'528.3/+20', '532.3/-20'};

% Plot the two observers with the default 590-nm test light
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors(1:2), l1, 'Case 1: no OD variation');
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors(1:2), l2, 'Case 2: 10% OD variation')
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments3,...
    colors(1:2), l3, 'Case 2: 20% OD variation')

% Plot L and M cone fundamentals for the two sets of observers
% The first number is the observer/lambda max, the
% second is the version of the observer (OD/no OD).
T_cones_11 = findConeFundamentals(lambdaMaxes(:,1),...
    dphotopigments1(:,1));
T_cones_12 = findConeFundamentals(lambdaMaxes(:,1),...
    dphotopigments2(:,1));
T_cones_13 = findConeFundamentals(lambdaMaxes(:,1),...
    dphotopigments3(:,1));
T_cones_21 = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments1(:,2));
T_cones_22 = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments2(:,2));
T_cones_23 = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments3(:,2));

figure; clf; hold on;
plot(wls,T_cones_11(1,:),'r','LineWidth',2);
plot(wls,T_cones_11(2,:),'r','LineWidth',2);
plot(wls,T_cones_21(1,:),'g','LineWidth',2);
plot(wls,T_cones_21(2,:),'g','LineWidth',2);
xlabel('Wavelength');
legend({'558.9', '528.3', '558.9', '532.3'});
title ('Compare cones - no OD variation');

figure; clf; hold on;
plot(wls,T_cones_12(1,:),'r','LineWidth',2);
plot(wls,T_cones_12(2,:),'r','LineWidth',2);
plot(wls,T_cones_22(1,:),'g','LineWidth',2);
plot(wls,T_cones_22(2,:),'g','LineWidth',2);
xlabel('Wavelength');
legend({'558.9/+10', '528.3/+10', '558.9/-10', '532.3/-10'});
title ('Compare cones - 10 percent OD variation');

figure; clf; hold on;
plot(wls,T_cones_13(1,:),'r','LineWidth',2);
plot(wls,T_cones_13(2,:),'r','LineWidth',2);
plot(wls,T_cones_23(1,:),'g','LineWidth',2);
plot(wls,T_cones_23(2,:),'g','LineWidth',2);
xlabel('Wavelength');
legend({'558.9/+20', '528.3/+20', '558.9/-20', '532.3/-20'});
title ('Compare cones - 20 percent OD variation');

% Find wavelength of maximum difference in each condition. Consider
% differences in both L and M cones above 540nm
% No OD variation
[d1, i1] = max(abs(T_cones_11(2,long_ind) - T_cones_21(2,long_ind)));
[d2, i2] = max(abs(T_cones_11(1,long_ind) - T_cones_21(1,long_ind)));
if d1 > d2
    ind1 = i1;
else
    ind1 = i2;
end
maxDiffWl1 = wls(ind1 + 160);

% 10% OD variation
[d1, i1] = max(abs(T_cones_12(2,long_ind) - T_cones_22(2,long_ind)));
[d2, i2] = max(abs(T_cones_12(1,long_ind) - T_cones_22(1,long_ind)));
if d1 > d2
    ind2 = i1;
else
    ind2 = i2;
end
maxDiffWl2 = wls(ind2 + 160);

% 20% OD variation
[d1, i1] = max(abs(T_cones_13(2,long_ind) - T_cones_23(2,long_ind)));
[d2, i2] = max(abs(T_cones_13(1,long_ind) - T_cones_23(1,long_ind)));
if d1 > d2
    ind3 = i1;
else
    ind3 = i2;
end
maxDiffWl3 = wls(ind3 + 160);

% Use the wavelengths of max difference as test wavelengths to separate
% the two observers
% No OD variation
plotTitle1 = sprintf('No Optical Density Variation, test light to %g nm',...
    maxDiffWl1);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors, l1, plotTitle1, 'test', maxDiffWl1);

% 10% OD variation
plotTitle2 = sprintf('10 percent Optical Density Variation, test light to %g nm',...
    maxDiffWl2);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors, l2, plotTitle2, 'test', maxDiffWl2);

% 20% OD variation
plotTitle3 = sprintf('20 percent Optical Density Variation, test light to %g nm',...
    maxDiffWl3);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments3,...
    colors, l3, plotTitle3, 'test', maxDiffWl3);

% Plot the overlapping observers (10% OD variation) with various test
% lights to see if any light can separate the two
test_wls = [540 560 580 600 620];
for i = 1:length(test_wls)
    theTitle = sprintf('Overlapping observers, test light = %g nm', test_wls(i)); 
    plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
        colors, l2, theTitle, 'test', test_wls(i));
end

%% Example 2 - L cone and OD variation
% lambdaMaxes +-2nm from average
lambdaMaxes = [
    [556.9 530.3 420.7]'...
    [560.9 530.3 420.7]'];
% Optical density +- 10%, 20% from average
dphotopigments1 = [...
    [0 0 0]'...
    [0 0 0]'];
dphotopigments2 = [...
    [15 15 0]'...
    [-15 -15 0]'];

% Legends
l1 = {'556.9', '560.9'};
l2 = {'556.9/+15', '560.9/-15'};

% Plot the two observers with the default 590-nm test light
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors(1:2), l1, 'Case 1: no OD variation');
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors(1:2), l2, 'Case 2: 15% OD variation')

% Plot L and M cone fundamentals for the two sets of observers
% The first number is the observer/lambda max, the
% second is the version of the observer (OD/no OD).
T_cones_11 = findConeFundamentals(lambdaMaxes(:,1),...
    dphotopigments1(:,1));
T_cones_12 = findConeFundamentals(lambdaMaxes(:,1),...
    dphotopigments2(:,1));
T_cones_21 = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments1(:,2));
T_cones_22 = findConeFundamentals(lambdaMaxes(:,2),...
    dphotopigments2(:,2));

figure; clf; hold on;
plot(wls,T_cones_11(1,:),'r','LineWidth',2);
plot(wls,T_cones_11(2,:),'r','LineWidth',2);
plot(wls,T_cones_21(1,:),'g','LineWidth',2);
plot(wls,T_cones_21(2,:),'g','LineWidth',2);
xlabel('Wavelength');
legend({'556.9', '530.3', '560.9', '530.3'});
title ('Compare cones - no OD variation');

figure; clf; hold on;
plot(wls,T_cones_12(1,:),'r','LineWidth',2);
plot(wls,T_cones_12(2,:),'r','LineWidth',2);
plot(wls,T_cones_22(1,:),'g','LineWidth',2);
plot(wls,T_cones_22(2,:),'g','LineWidth',2);
xlabel('Wavelength');
legend({'556.9/+15', '530.3/+15', '560.9/-15', '530.3/-15'});
title ('Compare cones - 15 percent OD variation');

% Find wavelength of maximum difference in each condition. Consider
% differences in both L and M cones above 540nm
% No OD variation
[d1, i1] = max(abs(T_cones_11(2,long_ind) - T_cones_21(2,long_ind)));
[d2, i2] = max(abs(T_cones_11(1,long_ind) - T_cones_21(1,long_ind)));
if d1 > d2
    ind1 = i1;
else
    ind1 = i2;
end
maxDiffWl1 = wls(ind1 + 160);

% 15% OD variation
[d1, i1] = max(abs(T_cones_12(2,long_ind) - T_cones_22(2,long_ind)));
[d2, i2] = max(abs(T_cones_12(1,long_ind) - T_cones_22(1,long_ind)));
if d1 > d2
    ind2 = i1;
else
    ind2 = i2;
end
maxDiffWl2 = wls(ind2 + 160);

% Use the wavelengths of max difference as test wavelengths to separate
% the two observers
% No OD variation
plotTitle1 = sprintf('No Optical Density Variation, test light to %g nm',...
    maxDiffWl1);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors, l1, plotTitle1, 'test', maxDiffWl1);

% 15% OD variation
plotTitle2 = sprintf('15 percent OD Variation, test light to %g nm',...
    maxDiffWl2);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors, l2, plotTitle2, 'test', maxDiffWl2);


% Plot the overlapping observers (15% OD variation) with various test
% lights to see if any light can separate the two
test_wls = [540 560 580 600 620];
for i = 1:length(test_wls)
    theTitle = sprintf('Overlapping observers, test light = %g nm', test_wls(i)); 
    plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
        colors, l2, theTitle, 'test', test_wls(i));
end