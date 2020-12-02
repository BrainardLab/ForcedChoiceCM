function limMatrix = qpTestWlParamRangeCheck(boundarySize,varargin)
% Finds approximate limits for matching range with a given test wavelength 
%
% Syntax:
%   qpTestWlParamRangeCheck(boundarySize)
%
% Description:
%    This function draws 100 observers and computes their nominal Rayleigh 
%    matches using a series of test wavelengths. The program calculates the
%    mean and standard deviation of matches for each test wavelength, then
%    uses these to define minimum and maximum limits for the matching
%    range. 
%
% Inputs:
%    boundarySize      -Numeric scalar specifying how many standard
%                       deviations away from the mean the limits should be
%                       set.
%
% Outputs:
%    limMatrix         -length(testWls) x 5 matrix for storing limits on 
%                       stimulus parameters. Each row represents a given  
%                       test wavelength and the limits which are associated 
%                       with it. The columns are arranged as follows: 
%                       [test wl, min lambda, max lambda, min test 
%                       intensity, max test intensity]. 
%
% Optional key-value pairs:
%    'nObservers'        -Number of simulated observers to test. Default is
%                         100.
%    'baseConeParams'    -Eight-element numeric vector of individual
%                         difference parameters: lens pigment density,
%                         macular pigment density, L/M/S photopigment
%                         densities, and L/M/S peak spectral sensitivities
%                         (lambda max). The values are used as means for
%                         observer sampling. Default is zeros(1,8).
%    'coneParamsToVary'  -Eight-element numeric vector of ones and zeros
%                         indicating which individual difference parameters
%                         should be varied (the noise parameter is excluded).
%                         Default is [0 0 1 1 0 1 1 0].
%    'p1Wl'              -Integer value of desired wavelength for the first
%                         primary. Default is 670.
%    'p2Wl'              -Integer value of desired wavelength for the second
%                         primary. Default is 560.
%    'testWls'           -Integer or numeric vector of desired wavelengths
%                         for the test light. Default is [570 590 610 630
%                         650].
%    'fieldSize'         -Integer field size, in degrees. Default is 2.
%    'age'               -Integer age for simulated observer. Default is
%                         32.
%    'opponentParams'    -4-element vector with opponent contrast
%                         parameters. (1) is the luminance weight, (2) is
%                         the RG weight, (3) is the BY weight, and (4) is
%                         the baseline noise standard deviation. Default is
%                         [40.3908  205.7353   62.9590    1.0000].
%    'p1Scale'           -Numerical scale factor for the first primary
%                         light, between 0 and 1. Default is 1.
%    'p2Scale'           -Numerical scale factor for the second primary
%                         light, between 0 and 1. Default is 0.02.
%    'testScale'         -Numerical scale factor for the test light,
%                         between 0 and 1. Default is 0.5.
%    'S'                 -Wavelength sampling for cone calculations, in the
%                         form [start increment numTerms]. Default is
%                         [380 2 201];
%    'makePlot'          -Logical. If true, makes a plot of the matches and
%                         limits.

% As of now, returns the following limit matrix for default settings:
% limMatrix = [570.0000,0.0090,0.1128,0.0410,0.0449;...
%              590.0000,0.0866,0.4278,0.0510,0.0673;...
%              610.0000,0.3140,0.7525,0.0809,0.1223;...
%              630.0000,0.6397,0.9215,0.1804,0.2523;...
%              650.0000,0.8814,0.9799,0.5329,0.6267];

% History
%    11/12/20   dce    -Wrote it
%    11/21/20   dce    -Changed to a function
%
% Example: qpTestWlParamRangeCheck(8)

% Parse input
p = inputParser;
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('nObservers',100,@(x)(isnumeric(x)));
p.addParameter('coneParamsToVary',[0 0 1 1 0 1 1 0],@(x)(isnumeric(x)));
p.addParameter('baseConeParams',zeros(1,8),@(x)(isnumeric(x)));
p.addParameter('p1Wl',670,@(x)(isnumeric(x)));
p.addParameter('p2Wl',560,@(x)(isnumeric(x)));
p.addParameter('testWls',[570 590 610 630 650],@(x)(isnumeric(x)));
p.addParameter('p1Scale',1,@(x)(isnumeric(x)));
p.addParameter('p2Scale',0.02,@(x)(isnumeric(x)));
p.addParameter('testScale',0.5,@(x)(isnumeric(x)));
p.addParameter('makePlot',true,@(x)(islogical(x)));
p.parse(varargin{:});

% Sample observers
observerParams = sampleRayleighObservers(p.Results.nObservers,...
    p.Results.baseConeParams,p.Results.coneParamsToVary);
observers = cell(1,p.Results.nObservers);
for i = 1:p.Results.nObservers
    observers{i} = genRayleighObserver('coneVec',observerParams(i,:),'age',...
        p.Results.age,'fieldSize',p.Results.fieldSize,'S',p.Results.S,...
        'opponentParams',p.Results.opponentParams);
end

% Generate spds
wls = SToWls(p.Results.S);
p1Spd = zeros(length(wls),1);
p1Spd(wls==p.Results.p1Wl) = 1*p.Results.p1Scale;

p2Spd = zeros(length(wls),1);
p2Spd(wls==p.Results.p2Wl) = 1*p.Results.p2Scale;

testSpds = zeros(length(wls),length(p.Results.testWls));
for i = 1:length(p.Results.testWls)
    testSpds(wls==p.Results.testWls(i),i) = 1*p.Results.testScale;
end

% Calculate nominal (analytic) matches for each observer at each test light.
% Store nominal matches in the third dimension of the array in the format
% [lambda testIntensity];
observerMatches = zeros(length(p.Results.testWls),p.Results.nObservers,2);
for i = 1:length(p.Results.testWls)
    for j = 1:p.Results.nObservers
        [~,~,testIntensity,lambda] = ...
            computePredictedRayleighMatch(p1Spd,p2Spd,testSpds(:,i),...
            observers{j},'addDarkSpd',false,'S',p.Results.S);
        observerMatches(i,j,1) = lambda;
        observerMatches(i,j,2) = testIntensity;
    end
end

% Calculate means and standard deviations
lambdaMeans = zeros(1,length(p.Results.testWls));
lambdaSDs = zeros(1,length(p.Results.testWls));
tIMeans = zeros(1,length(p.Results.testWls));
tISDs = zeros(1,length(p.Results.testWls));

for i = 1:length(p.Results.testWls)
    lambdas = observerMatches(i,:,1);
    lambdaMeans(i) = mean(lambdas);
    lambdaSDs(i) = std(lambdas);
    
    testIntensities = observerMatches(i,:,2);
    tIMeans(i) = mean(testIntensities);
    tISDs(i) = std(testIntensities);
end

% Get ranges based on a given number of standard deviations in any direction
lambdaMins = lambdaMeans-lambdaSDs*boundarySize;
lambdaMins(lambdaMins < 0) = 0;
lambdaMaxes = lambdaMeans+lambdaSDs*boundarySize;
tIMins = tIMeans-tISDs*boundarySize;
tIMins(tIMins < 0) = 0; 
tIMaxes = tIMeans+tISDs*boundarySize;
limMatrix = [p.Results.testWls',lambdaMins',lambdaMaxes',tIMins',tIMaxes'];

% Plot results
if p.Results.makePlot
    plotColors = 'rgbymck'; % Colors for plotting
    figure();
    hold on;
    for i = 1:length(p.Results.testWls)
        for j = 1:p.Results.nObservers
            % Cloud of points
            plot(observerMatches(i,j,1),observerMatches(i,j,2),'. ','Color',...
                plotColors(i));
        end
        % 3SD rectangle
        rectangle('Position',[lambdaMeans(i)-lambdaSDs(i)*3,tIMeans(i)-tISDs(i)*3,...
           lambdaSDs(i)*6,tISDs(i)*6],'EdgeColor','k');
        
        % Rectangle of computed limits
        rectangle('Position',[lambdaMeans(i)-lambdaSDs(i)*boundarySize...
            ,tIMeans(i)-tISDs(i)*boundarySize,lambdaSDs(i)*2*boundarySize,...
            tISDs(i)*2*boundarySize],'EdgeColor','r');
    end
    title('Distribution of Observer Matches');
    xlabel('Lambda');
    ylabel('Test Intensity');
end
end