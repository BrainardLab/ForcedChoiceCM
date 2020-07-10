function [primary_redder,test_brighter,isMatch,differenceVector] = ...
    observerRayleighDecision(observer,primarySpd,testSpd,varargin) 
% Simulated observer decision making for a Rayleigh match experiment 
%
% Syntax:
%   observerRayleighDecision(T_observer, primarySpd, testSpd); 
%
% Description:
%    Takes in a simulated observer's cone fundamentals and a pair of
%    spds (primary and test). 
%
%    Computes the opponent contrast of the two spectra based on the
%    observer's cone fundamentals. Using the contrast values of individual
%    channels, determines whether the test light should be made brighter or
%    dimmer and whether the mixing light should be shifted towards the red
%    primary or the green primary. Also determines whether the difference
%    between lights is below the matching threshold, which is a scalar 
%    multiple of the observer standard deviation. Noise can be added as a
%    key-value pair. 
%   
% Inputs:
%     observer          -Struct with observer parameters. Must contain the
%                        following fields: colorDiffParams, T_cones. 
%     primarySpd        -201x1 vector representation of the primary spd
%     testSpd           -201x1 vector represetation of the test spd  
%
% Outputs:
%     primary_redder    -Logical indicating whether primary ratio should 
%                        be shifted towards a redder value 
%     test_brighter     -Logical indicating whether test intensity should  
%                        be increased. 
%     isMatch           -Logical indicating whether the primary and test
%                        lights are below the matching threshold.
%
% Optional key-value pairs:
%     thresholdScale    -Number >0 which determines which scalar multiple 
%                        of the observer standard deviation should be used
%                        as the matching threshold. Default is 1. 
%     baseThreshold     -Number >0. When the observer is noiseless, 
%                        thresholdScale is multiplied by baseThreshold 
%                        instead of by the observer standard deviation. 
%                        Default is 0.02, an approximation of typical 
%                        observer standard deviation.

% History:
%   06/02/20  dce       Wrote initial code
%   06/03/20  dce       Added noise, style edits
%   06/09/20  dce       Added match threshold calculation  
%   06/11/20  dce       Added baseThreshold parameter
%   07/01/20  dce       Added check on observer's wavelength sampling 
%   07/06/20  dce       Got rid of "noisy" key-value pair, do noise based
%                       on SD only

% Parse input 
p = inputParser;
p.addParameter('thresholdScale',1, @(x)(isnumeric(x)));
p.addParameter('baseThreshold',0.02,@(x)(isnumeric(x)));
p.parse(varargin{:});

% Spds generated for OLRayleighMatch will always be sampled over S = 
% [380 2 201]. Check if the simulated observer follows this convention.
if ~all(observer.S == [380 2 201])
    error('Observer wavelength sampling does not follow OneLight convention'); 
end 

% Cone responses for the given spectra
test_LMS = observer.T_cones*testSpd; 
primary_LMS = observer.T_cones*primarySpd; 

% Opponent contrasts for the given spectra (primary mixing light relative 
% to test light). The result has the form [LUM; RG; BY].
opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
    test_LMS,primary_LMS);

% Add noise
% Sample three values from a Gaussian distribution, and add them to the
% opponentContrast vector
noiseSd = observer.colorDiffParams.noiseSd; % Gaussian standard deviation
if noiseSd ~= 0
    noise = normrnd(0,noiseSd,3,1); 
    opponentContrast = opponentContrast+noise; 
end 

% Check luminance 
if opponentContrast(1) > 0 % Primary mixture is brighter than test 
    test_brighter = true;
else 
    test_brighter = false; 
end 

% Check R/G ratio
if opponentContrast(2) < 0 % Primary mixture is greener than test
    primary_redder = true; 
else 
    primary_redder = false; 
end

% Define matching threshold 
if noiseSd ~= 0
    threshold = noiseSd*p.Results.thresholdScale;  
else 
    threshold = p.Results.baseThreshold*p.Results.thresholdScale;
end

% Check whether we're below the matching threshold  
differenceVector = norm(opponentContrast(1),opponentContrast(2)); 
if differenceVector < threshold
    isMatch = true; 
else
    isMatch = false; 
end 
end 