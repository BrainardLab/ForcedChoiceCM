function [primary_redder,test_brighter,isMatch,differenceVector] = ...
    observerRayleighDecision(observer,primarySpd,testSpd,varargin) 
% Simulated observer decision making for a Rayleigh match experiment 
%
% Syntax:
%   observerRayleighDecision(observer, primarySpd, testSpd); 
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
%     primarySpd        -nWlsx1 vector representation of the primary spd
%     testSpd           -nWlsx1 vector represetation of the test spd  
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
%     'noiseScale'      -Number >=0 which determines which scalar multiple
%                        of the opponent noise Sd should be used as the
%                        observer noise Sd. Default is 1.
%     'thresholdScale'  -Number >0 which determines which scalar multiple 
%                        of the observer noise should be used as the 
%                        matching threshold. Default is 0.5.
%     'refSpd'          -nWlsx1 vector represetation of the reference spd.
%                        If this is nonempty, the program computes the 
%                        opponent contrast of both primary mixture and test 
%                        light relative to the reference, rather than their 
%                        opponent contrast relative to one another. Default 
%                        is [].
% History:
%   06/02/20  dce       Wrote initial code
%   06/03/20  dce       Added noise, style edits
%   06/09/20  dce       Added match threshold calculation  
%   06/11/20  dce       Added baseThreshold parameter
%   07/01/20  dce       Added check on observer's wavelength sampling 
%   07/06/20  dce       Got rid of "noisy" key-value pair, do noise based
%                       on SD only
%   08/05/20  dce       Distinguished opponent noise and observer noise.
%   11/15/20  dce       Added option to compute opponent contrast relative 
%                       to a reference, changed so opponent contrast of
%                       test is computed relative to primary mixture.

% Parse input 
p = inputParser;
p.addParameter('noiseScale',1, @(x)(isnumeric(x)));
p.addParameter('thresholdScale',0.5, @(x)(isnumeric(x)));
p.addParameter('refSpd',[], @(x)(isnumeric(x)));
p.parse(varargin{:});

% Cone responses for the given spectra
test_LMS = observer.T_cones*testSpd; 
primary_LMS = observer.T_cones*primarySpd; 

% Opponent contrasts for the given spectra. The result has the form [LUM; RG; BY].
if isempty(p.Results.refSpd)
    % No reference light provided - return opponent contrast of test
    % relative to primary mixture.
    opponentContrastDiff = LMSToOpponentContrast(observer.colorDiffParams,...
        primary_LMS,test_LMS);
else
    % Reference light provided - return opponent contrast of test -
    % opponent contrast of primary mixture.
    ref_LMS = observer.T_cones*p.Results.refSpd;
    pOpponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
        ref_LMS,primary_LMS);
    tOpponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
        ref_LMS,test_LMS);
    opponentContrastDiff = tOpponentContrast - pOpponentContrast;
end 

% Add noise to opponent contrast by sampling three values from a Gaussian
% distribution.
opponentNoiseSd = observer.colorDiffParams.noiseSd; 
observerNoiseSd = p.Results.noiseScale*opponentNoiseSd;
if observerNoiseSd ~= 0 % The simulated matching is noisy
    noise = normrnd(0,observerNoiseSd,3,1); 
    opponentContrastDiff = opponentContrastDiff+noise; 
end 

% Check luminance. If negative, the test light needs to be made brighter. 
if opponentContrastDiff(1) < 0 % Test is less bright than primary mixture 
    test_brighter = true;
else 
    test_brighter = false; 
end 

% Check R/G ratio. If negative, the primary mixture needs to be made
% redder.
if opponentContrastDiff(2) < 0 % Test is redder than primary mixture 
    primary_redder = false; 
else 
    primary_redder = true; 
end

% Define matching threshold 
if observerNoiseSd ~= 0
    threshold = observerNoiseSd*p.Results.thresholdScale;  
else 
    threshold = opponentNoiseSd*p.Results.thresholdScale;
end

% Check whether we're below the matching threshold  
differenceVector = norm(opponentContrastDiff(1:2)); 
if differenceVector < threshold
    isMatch = true; 
else
    isMatch = false; 
end 
end 