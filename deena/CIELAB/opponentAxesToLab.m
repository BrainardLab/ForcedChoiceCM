% This script is designed to find the appropriate scalars for opponent
% contrast metrics to bring the opponent contrast space in line with 
% CIELAB scaling. This would mean that an equal step size along either of  
% the three axes (luminance, rg, and by) corresponds to the same perceptual
% difference.
%
% The script starts with a set of reference LMS coordinates, and it samples
% a set of cooordinates at equal CIELAB distance away from the reference
% coordinates (using a sphere). The reference and sampled coordinates are
% used in a parameter search function, which searches for the optimal
% luminance, rg, and by opponenet channel scalars. The goal is to find
% scalars that lead the opponent contrast of two coordinates to be equal in
% magnitude to the CIELAB distance between them. In this script, both
% distances are set equal to the observer standard deviation. 

% History:
%   07/29/20  dce, dhb   Wrote it.

%% Generate a set of LMS coordinates that are at a uniform distance from 
%% the reference in LAB space
% Initial parameters 
refLMS = [1 1 1]';  % Reference LMS coordinates (set later)
noiseSD = 1;        % Observer noise standard deviation (step size/radius)

% Convert reference to Lab coordinates, using itself as white point
[refLab,~,refXYZ] = LMSToLAB(refLMS,refLMS); 
if ~all(refLab==[100 0 0]')
    error('RefLab not computed properly');
end 

% Sample Lab points on a sphere around refLab. 
nTheta = 10;  % Number of points to sample around the azimuthal theta 
nPhi = 10;    % Number of points to sample around the elevation direction
nPointsOnSphere = nTheta*nPhi; 

% Generate a sphere, and adjust so it is centered at the refLab
% coordinates. NoiseSD is used as the radius. 
sphereLab = SphereGenerate(nTheta,nPhi,1);
sphereLab(1,:) = sphereLab(1,:)+refLab(1); 
sphereLab(2,:) = sphereLab(2,:)+refLab(2); 
sphereLab(3,:) = sphereLab(3,:)+refLab(3); 

% Convert sphereLab points back to LMS
sphereLMS = zeros(size(sphereLab));
for i = 1:nPointsOnSphere
    sphereLMS(:,i) = LABToLMS(sphereLab(:,i),refXYZ);
end 

%% Parameter search for opponent metric scalings 
% We want to choose our opponent metric scalings so that the opponent
% space distance between refLMS and sphereLMS is as close to the noise SD 
% as we can make it.

% Fit three params: luminance weight (1), rg weight (1), and by weight(1)
initialParams = noiseSD*[1 1 1]/0.02;  % Initial parameters  

% fmincon options  
lb = -Inf*ones(3,1);     % Lower bounds
ub = Inf*ones(3,1);      % Upper bounds
options = optimset('fmincon');
options = optimset(options,'Diagnostics','off','Display','iter',...
    'LargeScale','off','Algorithm','active-set');

% Run the test 
bestODParams = fmincon(@(x)opponentParamsErr(x,refLMS,sphereLMS,noiseSD),...
    initialParams,[],[],[],[],lb,ub,[],options);
err = opponentParamsErr(bestODParams,refLMS,sphereLMS,noiseSD);