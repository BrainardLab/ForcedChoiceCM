function [opponentContrastParams,error,opponentSphere] = opponentAxesToLab(noiseSD)
% Syntax:
%   opponentParamsErr(opponentParams,refLMS,shiftedLMS,noiseSD)
%
% Description:
%    Finds the optimal opponent contrast channel weights (luminance
%    weighting, rg weighting, and by weighting) that bring the opponent
%    space scaling in line with CIELAB space. The program samples a sphere
%    of points that are a CIELAB distance of 1 from a reference, then finds
%    opponent weighting parameters that lead the opponent contrast to be
%    equal to the selected noise standard deviation on average. The search
%    is conducted using fmincon.
%
%    With the current settings, the program returns ideal params of
%    [0.8078 4.1146 1.2592 0.0200], which are used as defaults
%    throughout the simulation programs.
%
% Inputs:
%    noiseSD        -Numerical measure of noise standard deviation
%
% Outputs:
%    opponentContrastParams -4-element vector containing optimal opponent
%                            contrast parameters. (1) is luminance weight,
%                            (2) is rg weight, (3) is by weight, and (4) is
%                            observer noise.
%    error                  -rms error of the difference between opponent
%                            contrast norm and the noise standard deviation
%
% Optional key-value pairs:
%    None

% History:
%   07/29/20  dce, dhb   Wrote it.
%   08/03/20  dce        Made a function, added plotting
%   08/10/20  dce        Added 3D plots, passed in cones to LMS-LAB
%                        conversion functions

%% Generate a set of LMS coordinates that are at a uniform distance from
%% the reference in LAB space
% Check input
if noiseSD <= 0
    error('Noise must be greater than 0');
end

% Load data files
coneData = load('T_cones_ss2');
xyzData = load('T_xyz1931');
spdData = load('spd_D65');

% Calculate reference LMS as the response to D65 illuminant
spd_D65 = SplineSpd(spdData.S_D65,spdData.spd_D65,coneData.S_cones_ss2);
refLMS = coneData.T_cones_ss2*spd_D65;
refLMS = refLMS/max(refLMS(:));

% Convert reference to Lab coordinates, using itself as white point
[refLab,~,refXYZ] = LMSToLAB(refLMS,refLMS,'T_cones',coneData.T_cones_ss2,...
    'S_cones',coneData.S_cones_ss2,'T_xyz',xyzData.T_xyz1931,'S_xyz',...
    xyzData.S_xyz1931);
if ~all(refLab==[100 0 0]')
    error('RefLab not computed properly');
end

% Sample Lab points on a sphere around refLab.
nTheta = 10;  % Number of points to sample around the azimuthal theta
nPhi = 10;    % Number of points to sample around the elevation direction
nPointsOnSphere = nTheta*nPhi;

% Generate a unit sphere, and adjust so it is centered at the refLab
% coordinates. The radius of 1 corresponds to a step of 1 in the CIELAB
% space, which signifies a barely-detectable perceptual difference.
baseSphere = UnitSphereGenerate(nTheta,nPhi);
sphereLab = zeros(size(baseSphere));
sphereLab(1,:) = baseSphere(1,:)+refLab(1);
sphereLab(2,:) = baseSphere(2,:)+refLab(2);
sphereLab(3,:) = baseSphere(3,:)+refLab(3);

% Convert sphereLab points back to LMS
sphereLMS = zeros(size(sphereLab));
for i = 1:nPointsOnSphere
    sphereLMS(:,i) = LABToLMS(sphereLab(:,i),refXYZ,'T_cones',...
        coneData.T_cones_ss2,'S_cones',coneData.S_cones_ss2,'T_xyz',...
        xyzData.T_xyz1931,'S_xyz',xyzData.S_xyz1931);
end

%% Parameter search for opponent metric scalings
% We want to choose our opponent metric scalings so that the opponent
% space distance between refLMS and sphereLMS is as close to the noise SD
% as we can make it.

% Fit three params: luminance weight (1), rg weight (1), and by weight(1).
% The initial values are scaled to a rough order of magnitude.
initialParams = [1 1 1]*(noiseSD/0.02);

% fmincon options. Can add more bounds as desired
lb = -Inf*ones(3,1);     % Lower bounds
ub = Inf*ones(3,1);      % Upper bounds
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter',...
    'LargeScale','off','Algorithm','active-set');

% Run the test
bestODParams = fmincon(@(x)opponentParamsErr(x,refLMS,sphereLMS,noiseSD),...
    initialParams,[],[],[],[],lb,ub,[],options);
opponentContrastParams = [bestODParams noiseSD];
errScalar = 100;
error = opponentParamsErr(bestODParams,refLMS,sphereLMS,noiseSD)/errScalar;

% Plot to see how close to a sphere the opponent points are
% Generate params struct
colorDiffParams = getColorDiffParams([bestODParams noiseSD]);

% Find opponent contrast vector
opponentContrasts = zeros(1,nPointsOnSphere);
for i = 1:nPointsOnSphere
    opponentContrast = LMSToOpponentContrast(colorDiffParams,refLMS,...
        sphereLMS(:,i));
    opponentContrasts(i) = norm(opponentContrast);
end

% Find opponent sphere coordinates in 3D
opponentCenter = LABToLMS(refLMS,refXYZ);  % Center of opponent sphere
opponentSphere(1,:) = baseSphere(1,:)+opponentCenter(1);
opponentSphere(2,:) = baseSphere(2,:)+opponentCenter(2);
opponentSphere(3,:) = baseSphere(3,:)+opponentCenter(3);
for i = 1:nPointsOnSphere  % Scale points by opponent radius
    opponentSphere(:,i) = opponentSphere(:,i)*opponentContrasts(i);
end
plotResults = false;
if plotResults
    % 2D plot
    figure();
    plot(1:nPointsOnSphere,opponentContrasts,'ro ');
    ylim([0 2*noiseSD]);
    refline(0,noiseSD);
    legend('Opponent Contrasts','Noise Standard Deviation');
    title(['Opponent Contrast Vector Lengths, Noise = ' num2str(noiseSD)]);
    xlabel('Point');
    ylabel('Vector Length');
    
    % 3D plots
    figure();
    plot3(sphereLab(1,:),sphereLab(2,:),sphereLab(3,:),'r* ');
    title(['CIELAB Sphere, Noise = ' num2str(noiseSD)]);
    xlim([refLab(1)-1.25,refLab(1)+1.25]);
    ylim([refLab(2)-1.25,refLab(2)+1.25]);
    zlim([refLab(3)-1.25,refLab(3)+1.25]);
    xlabel('x coords');
    ylabel('y coords');
    zlabel('z coords');
    
    figure();
    plot3(opponentSphere(1,:),opponentSphere(2,:),opponentSphere(3,:),'r* ');
    title(['Opponent Sphere, Noise = ' num2str(noiseSD)]);
    xlim([opponentCenter(1)-1.25*noiseSD,opponentCenter(1)+1.25*noiseSD]);
    ylim([opponentCenter(2)-1.25*noiseSD,opponentCenter(2)+1.25*noiseSD]);
    zlim([opponentCenter(3)-1.25*noiseSD,opponentCenter(3)+1.25*noiseSD]);
    xlabel('x coords');
    ylabel('y coords');
    zlabel('z coords');
end
end