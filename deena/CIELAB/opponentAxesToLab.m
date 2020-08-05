function [opponentContrastParams,error] = opponentAxesToLab(noiseSD)
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

%% Generate a set of LMS coordinates that are at a uniform distance from
%% the reference in LAB space
% Check input
if noiseSD <= 0 
    error('Noise must be greater than 0');
end 

% Initial parameters
refLMS = [1 1 1]';  % Reference LMS coordinates (set later)

% Convert reference to Lab coordinates, using itself as white point
[refLab,~,refXYZ] = LMSToLAB(refLMS,refLMS);
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
plotResults = false;
if plotResults
    % Generate output struct
    colorDiffParams = struct();
    colorDiffParams.type = 'opponentContrast';
    colorDiffParams.LMRatio = 2;
    colorDiffParams.lumWeight = bestODParams(1);
    colorDiffParams.rgWeight = bestODParams(2);
    colorDiffParams.byWeight = bestODParams(3);
    colorDiffParams.M = GetOpponentContrastMatrix(colorDiffParams);
    opponentContrasts = zeros(1,nPointsOnSphere);
    for i = 1:nPointsOnSphere
        opponentContrast = LMSToOpponentContrast(colorDiffParams,refLMS,...
            sphereLMS(:,i));
        opponentContrasts(i) = norm(opponentContrast);
    end
    figure();
    plot(1:nPointsOnSphere,opponentContrasts,'ro ');
    ylim([0 2*noiseSD]);
    refline(0,noiseSD);
    legend('Opponent Contrasts','Noise Standard Deviation');
    title(['Opponent Contrast Vector Lengths, Noise = ' num2str(noiseSD)]);
    xlabel('Point');
    ylabel('Vector Length');
end
end