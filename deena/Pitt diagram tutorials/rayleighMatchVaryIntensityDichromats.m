% Script for investigating how changing different lambda max/optical
% density parameters for dichromats affects simulated Rayleigh matches.
% Plots the ranges of these parameters across lambda max/optical density
% variation, then implements two case studies where optical density 
% variations obscure lambda max variations (one for deuteranopes and one
% for protanopes). Identifies the wavelength with the largest difference in
% cone activations and uses this as a test light to separate the two cones.

%% Physiological ranges of parameters
% Cone lambda maxes for protanopes (M cones: 530.3+-6nm)
lambdaMaxesP = [  ...
    [524.3 524.3 420.7]' ...
    [526.3 526.3 420.7]'...
    [528.3 528.3 420.7]'...
    [530.3 530.3 420.7]'...
    [532.3 532.3 420.7]'...
    [534.3 534.3 420.7]'...
    [536.3 536.3 420.7]'];...
    
lambdaMaxesP_small = [  ...
    [524.3 524.3 420.7]' ...
    [530.3 530.3 420.7]'...
    [536.3 536.3 420.7]'];...
    
% Cone lambda maxes for deuteranopes (L cones: 558.9+-6nm)
lambdaMaxesD = [  ...
    [552.9 552.9 420.7]' ...
    [554.9 554.9 420.7]'...
    [556.9 556.9 420.7]'...
    [558.9 558.9 420.7]'...
    [560.9 560.9 420.7]'...
    [562.9 562.9 420.7]'...
    [564.9 564.9 420.7]'];

lambdaMaxesD_small = [  ...
    [552.9 552.9 420.7]' ...
    [558.9 558.9 420.7]'...
    [564.9 564.9 420.7]'];

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
[~, nump] = size(lambdaMaxesP);
[~, numd] = size(lambdaMaxesD);
colors = 'bgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmykbgrcmyk';

%% Case 1: Vary only lambda max
% Protanopes
for i = 1:numODs
    dphotopigments = repmat(ODs(:,i), 1, numd);
    theTitle = sprintf('Protanopes, OD = %g percent', ODs(1,i));
    theLegend = {'524.3', '526.3', '528.3', '530.3', '532.3', '534.3', '536.3'};
    plotRayleighMatches(lambdaMaxesP, lambdaMaxesP, dphotopigments,...
        colors(1:nump), theLegend, theTitle);
end

% Deuteranopes
for i = 1:numODs
    dphotopigments = repmat(ODs(:,i), 1, numd);
    theTitle = sprintf('Deuteranopes, OD = %g percent', ODs(1,i));
    theLegend = {'552.9', '554.9', '556.9', '558.9', '560.9', '562.9', '564.9'};
    plotRayleighMatches(lambdaMaxesD, lambdaMaxesD, dphotopigments,...
        colors(1:numd), theLegend, theTitle);
end

%% Case 2: Vary only OD
% Protanopes
for i = 1:nump
    lambdaMaxes = repmat(lambdaMaxesP(:,i), 1, numODs);
    theTitle = sprintf('OD variation, lambdaMax = %g nm', lambdaMaxesP(1,i));
    theLegend = {'-30', '-20', '-10', '0', '10', '20', '30'};
    plotRayleighMatches(lambdaMaxes, lambdaMaxes, ODs,...
        colors(1:numODs), theLegend, theTitle);
end

% Deuteranopes
for i = 1:numd
    lambdaMaxes = repmat(lambdaMaxesD(:,i), 1, numODs);
    theTitle = sprintf('OD variation, lambdaMax = %g nm', lambdaMaxesD(1,i));
    theLegend = {'-30', '-20', '-10', '0', '10', '20', '30'};
    plotRayleighMatches(lambdaMaxes, lambdaMaxes, ODs,...
        colors(1:numODs), theLegend, theTitle);
end

%% Case 3: Vary simultaneously
% Protanopes
lambdaMaxes = [];
dphotopigments = [];
for i = 1:length(lambdaMaxesP_small)
    lambdaMaxes = [lambdaMaxes, repmat(lambdaMaxesP_small(:,i), 1, length(ODs_small))];
    dphotopigments = [dphotopigments, ODs_small];
end
theLegend = {'524.3/-30', '524.3/-15', '524.3/0', '524.3/15', '524.3/30'...
    '530.3/-30', '530.3/-15', '530.3/0', '530.3/15', '530.3/30'...
    '536.3/-30', '536.3/-15', '536.3/0', '536.3/15', '536.3/30'}; 
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments,...
    colors, theLegend, 'Simultaneous Variation - Deuteranopes');

% Deuteranopes
lambdaMaxes = [];
dphotopigments = [];
for i = 1:length(lambdaMaxesD_small)
    lambdaMaxes = [lambdaMaxes, repmat(lambdaMaxesD_small(:,i), 1, length(ODs_small))];
    dphotopigments = [dphotopigments, ODs_small];
end
theLegend = {'552.9/-30', '552.9/-15', '552.9/0', '552.9/15', '552.9/30'...
    '558.9/-30', '558.9/-15', '558.9/0', '558.9/15', '558.9/30'...
    '564.9/-30', '564.9/-15', '564.9/0', '564.9/15', '564.9/30'}; 
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments,...
    colors, theLegend, 'Simultaneous Variation - Protanopes');

%% Case studies
% Show that two distinguishable lambdaMaxes can be made to overlap through
% OD variation, and try to separate them with test lights 

% Wavelength sampling
S = [380 1 401];
wls = SToWls(S);
long_ind = 161:length(wls); % Indices of wavelengths above 540nm

%% 1 - Protanopes 
% lambdaMaxes +-2nm from average
lambdaMaxes = [
    [528.3 528.3 420.7]'...
    [532.3 532.3 420.7]']; 
% Optical density +- 20% from average
dphotopigments1 = [...
    [0 0 0]'...
    [0 0 0]']; 
dphotopigments2 = [...
    [20 20 0]'...
    [-20 -20 0]']; 
l1 = {'528.3', '532.3'}; 
l2 = {'528.3/+20', '532.3/-20'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors(1:2), l1, 'Case 1: no OD variation');
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors(1:2), l2, 'Case 2: OD variation')

% Plot cone fundamentals for the two sets of observers
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
plot(wls,T_cones_21(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (l1);
title ('Compare M cones - no OD variation');

figure; clf; hold on;
plot(wls,T_cones_12(1,:),'r','LineWidth',2);
plot(wls,T_cones_22(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (l2);
title ('Compare M cones - OD variation');

% Find wavelength of maximum difference in each condition
[~, ind1] = max(abs(T_cones_11(1,long_ind) - T_cones_21(1,long_ind)));
maxDiffWl1 = wls(ind1 + 160);

[~, ind2] = max(abs(T_cones_12(1,long_ind) - T_cones_22(1,long_ind)));
maxDiffWl2 = wls(ind2 + 160);

% Use the wavelengths of max difference as test wavelengths to separate
% the two observers
plotTitle1 = sprintf('No Optical Density Variation, test light to %g nm',...
    maxDiffWl1);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors, l1, plotTitle1, 'test', maxDiffWl1);

plotTitle2 = sprintf('Optical Density Variation, test light to %g nm',...
    maxDiffWl2);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors, l2, plotTitle2, 'test', maxDiffWl2);

%% 2 - Deuteranopes 
% lambdaMaxes +-2nm from average
lambdaMaxes = [
    [557.9 557.9 420.7]'...
    [559.9 559.9 420.7]']; 
% Optical density +- 20% from average
dphotopigments1 = [...
    [0 0 0]'...
    [0 0 0]']; 
dphotopigments2 = [...
    [20 20 0]'...
    [-20 -20 0]']; 
l1 = {'557.9', '559.9'}; 
l2 = {'557.9/+20', '559.9/-20'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors(1:2), l1, 'Case 1: no OD variation');
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors(1:2), l2, 'Case 2: OD variation')

% Find cone fundamentals. The first number is the observer/lambda max, the
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
plot(wls,T_cones_21(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (l1);
title ('Compare L cones - no OD variation');

figure; clf; hold on;
plot(wls,T_cones_12(1,:),'r','LineWidth',2);
plot(wls,T_cones_22(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (l2);
title ('Compare L cones - OD variation');

% Find wavelength of maximum difference in each condition
[~, ind1] = max(abs(T_cones_11(1,long_ind) - T_cones_21(1,long_ind)));
maxDiffWl1 = wls(ind1 + 160);

[~, ind2] = max(abs(T_cones_12(1,long_ind) - T_cones_22(1,long_ind)));
maxDiffWl2 = wls(ind2 + 160);

% Use the wavelengths of max difference as test wavelengths to separate
% the two observers
plotTitle1 = sprintf('No Optical Density Variation, test light to %g nm',...
    maxDiffWl1);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors, l1, plotTitle1, 'test', maxDiffWl1);

plotTitle2 = sprintf('Optical Density Variation, test light to %g nm',...
    maxDiffWl2);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors, l2, plotTitle2, 'test', maxDiffWl2);

%% 3 - Deuteranopes - 2nm difference 
% lambdaMaxes +-1nm from average
lambdaMaxes = [
    [557.9 557.9 420.7]'...
    [559.9 559.9 420.7]']; 
% Optical density +- 20% from average
dphotopigments1 = [...
    [0 0 0]'...
    [0 0 0]']; 
dphotopigments2 = [...
    [20 20 0]'...
    [-20 -20 0]']; 
l1 = {'557.9', '559.9'}; 
l2 = {'557.9/+20', '559.9/-20'};
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors(1:2), l1, 'Case 1: no OD variation');
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors(1:2), l2, 'Case 2: OD variation')

% Find cone fundamentals. The first number is the observer/lambda max, the
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
plot(wls,T_cones_21(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (l1);
title ('Compare L cones - no OD variation');

figure; clf; hold on;
plot(wls,T_cones_12(1,:),'r','LineWidth',2);
plot(wls,T_cones_22(1,:),'g','LineWidth',2);
xlabel('Wavelength');
legend (l2);
title ('Compare L cones - OD variation');

% Find wavelength of maximum difference in each condition
[~, ind1] = max(abs(T_cones_11(1,long_ind) - T_cones_21(1,long_ind)));
maxDiffWl1 = wls(ind1 + 160);

[~, ind2] = max(abs(T_cones_12(1,long_ind) - T_cones_22(1,long_ind)));
maxDiffWl2 = wls(ind2 + 160);

% Use the wavelengths of max difference as test wavelengths to separate
% the two observers
plotTitle1 = sprintf('No Optical Density Variation, test light to %g nm',...
    maxDiffWl1);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments1,...
    colors, l1, plotTitle1, 'test', maxDiffWl1);

plotTitle2 = sprintf('Optical Density Variation, test light to %g nm',...
    maxDiffWl2);
plotRayleighMatches(lambdaMaxes, lambdaMaxes, dphotopigments2,...
    colors, l2, plotTitle2, 'test', maxDiffWl2);

