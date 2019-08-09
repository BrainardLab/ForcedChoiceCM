% Illustrate and test likelihood computation routines
%
% Description:
%

% History:
%   08/08/19  dhb  Wrote it.

%% Clear
clear; close all;

%% Set up parameters
S = [400 1 301];
params = ObserverVecToParams(zeros(8,1));
fieldSizeDegrees = 10;
ageInYears = 32;
pupilDiameterMM = 3;

%% Get cone fundamentals, just to make sure they look sensible
[~,T] = ComputeChoiceLikelihood(params,S,fieldSizeDegrees,ageInYears,pupilDiameterMM);
coneFunadmentalFig = figure; clf; hold on;
plot(SToWls(S),T(1,:)','r','LineWidth',2);
plot(SToWls(S),T(2,:)','g','LineWidth',2);
plot(SToWls(S),T(3,:)','b','LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Fundamental');


