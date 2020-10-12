function [predictedProportions] = qpPFRM(stimParamsVec,coneParamsVec,opponentParamsVec,noiseSD,S,p1Spd,p2Spd,tSpds,tWls,varargin)
% Psychometric function for Rayleigh matching
%
% Usage:
%     [predictedProportions] =
%     qpPFRM(stimParamsVec,coneParamsVec,opponentParamsVec,noiseSD,S,p1Spd,p2Spd,tSpds,tWls)
%
% Description:
%     Psychometric function for Rayleigh matching of luminance or 
%     chromaticity. Takes in a matrix of stimulus parameters which are used 
%     to generate a primary mixture and test light. Calculates the opponent  
%     contrast representation of the two lights for a specified observer, 
%     then returns predicted proportions for observer comparisons of 
%     either luminance or chromaticity.
%
% Input:
%     stimParamsVec      Matrix, with each row being a vector of stimulus
%                        parameters [lambda,test intensity,test peak wavelength].
%                        (lambda is the proportion of p1 in the primary
%                        mixture). 
%     coneParamsVec      Row vector of cone parameters for the observer.
%                        Consists of eight individual difference parameters 
%                        that follow the Asano model (2016). See  
%                        ObserverVecToParams for a full description. 
%     opponentParamsVec  4-element vector with obsever's opponent contrast 
%                        parameters. (1) is the luminance weight, (2) is 
%                        the rg weight,(3) is the by weight, and (4) is the 
%                        noise standard deviation.
%     noiseSD            Observer noise standard deviation.
%     S                  Wavelength sampling.
%     p1Spd              Spd for first primary in mixture (red), entered as 
%                        a column vector.
%     p2Spd              Spd for second primary in mixture (green), entered
%                        as a column vector.
%     tSpds              Matrix of spds for possible test lights at 
%                        different peak wavelengths. Has dimensions of spdLength 
%                        x nWls. 
%     tWls               Vector specifying which peak wavelengths are included
%                        in tSpds. Its length must equal the width of tSpds.
%
% Output:
%     predictedProportions   Matrix, where each row is a vector of 
%                            predicted proportions for each outcome.
%                            The first entry of each row is the "no" 
%                            response (test is less red/less bright than 
%                            primary), and the second is the "yes" response 
%                            (test is redder/brighter than primary). Note
%                            that this conventions differs from that of
%                            observerRayleighDecision.
%
% Optional key/value pairs
%     'judgeLum'             Compute proportions for a luminance judgement,  
%                            not a red/green judgement. Default is
%                            false.
%
% 10/08/20  dce  Wrote it 

%% Parse input
p = inputParser;
p.addRequired('stimParamsVec',@isnumeric);
p.addRequired('coneParamsVec',@isnumeric);
p.addRequired('opponentParamsVec',@isnumeric);
p.addRequired('noiseSD',@isnumeric);
p.addRequired('S',@isnumeric);
p.addRequired('p1Spd',@isnumeric);
p.addRequired('p2Spd',@isnumeric);
p.addRequired('tSpds',@isnumeric);
p.addRequired('tWls',@isnumeric);
p.addParameter('judgeLum',false,@(x)(islogical(x)));
p.parse(stimParamsVec,coneParamsVec,opponentParamsVec,noiseSD,S,p1Spd,...
    p2Spd,tSpds,tWls,varargin{:});

%% Input checks 
if length(tWls) ~= size(tSpds,2)
    error('Passed test wavelengths do not match passed test spds'); 
elseif ~all(ismember(stimParamsVec(:,3),tWls))
    error('Chosen test wavelengths not found in passed array');
elseif ~all(stimParamsVec(:,1) >= 0) || ~all(stimParamsVec(:,1) <= 1)
    error('Chosen lambda values must be between 0 and 1');
elseif ~all(stimParamsVec(:,2) >= 0) || ~all(stimParamsVec(:,2) <= 1)
    error('Chosen test intensity values must be between 0 and 1');
end 

%% Generate observer based on passed parameters 
observer = genRayleighObserver('coneVec',coneParamsVec,'opponentParams',...
    opponentParamsVec,'S',S); 

%% Loop over stimuli and get the proportions for each 
nStim = size(stimParamsVec,1);
predictedProportions = zeros(nStim,2);
for ii = 1:nStim
    % Extract the parameters
    lambda = stimParamsVec(ii,1);
    tI = stimParamsVec(ii,2); 
    testWl = stimParamsVec(ii,3); 
    
    % Generate the stimuli - primary mixture and test light 
    primary = lambda*p1Spd + (1-lambda)*p2Spd;
    test = tI*tSpds(:,tWls==testWl); 
    
    % Compute LMS responses of observer to the two lights
    pLMS = observer.T_cones * primary;
    tLMS = observer.T_cones * test;
    
    % Compute opponent contrast of test light relative to the primary
    % mixture (the primary is displayed for longer by default in the 
    % OneLight version of the experiment)
    opponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
        pLMS,tLMS);
    
    % Which channel are we looking at? Use this to set the mean of the
    % cumulative distribution.
    if p.Results.judgeLum
        mu = opponentContrast(1); % Luminance channel 
    else 
        mu = opponentContrast(2); % RG channel 
    end 
    
    % Generate probability distribution for the test light being
    % less bright/less red than primary 
    predictedProportions(ii,1) = normcdf(0,mu,noiseSD);
    
    % Fill in complement - probability of the test light being
    % brighter/redder than primary
    predictedProportions(ii,2) = 1-predictedProportions(ii,1);  
end

%% Don't allow complete certainty
predictedProportions(predictedProportions > 0.9999) = 0.9999;
predictedProportions(predictedProportions < 0.0001) = 0.0001;