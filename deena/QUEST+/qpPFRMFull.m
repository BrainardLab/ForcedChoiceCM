function [predictedProportions] = qpPFRMFull(stimParamsVec,coneParamsVec,...
    opponentParamsVec,noiseSD,S,p1Spd,p2Spd,tSpds,tWls,varargin)
% Psychometric function for Rayleigh matching
%
% Usage:
%     [predictedProportions] =
%     qpPFRMFull(stimParamsVec,coneParamsVec,opponentParamsVec,noiseSD,S,p1Spd,p2Spd,tSpds,tWls,lambdaRef)
%
%
% Description:
%     Psychometric function for Rayleigh matching of both luminance and  
%     red/green. Takes in a matrix of stimulus parameters which are used 
%     to generate a primary mixture and test light from provided initial 
%     spectra. Calculates the opponent contrast representation of the two 
%     lights for a specified observer, then returns predicted proportions 
%     for each of four possible outcomes: too red and too bright, too red 
%     and not too bright, too green and too bright, too green and not too 
%     bright. 
%
% Input:
%     stimParamsVec      Matrix, with each row being a vector of stimulus
%                        parameters [lambda,test intensity,test peak 
%                        wavelength].(lambda is the proportion of p1 in the 
%                        primary mixture). 
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
%                        different peak wavelengths. Has dimensions of  
%                        spdLength x nWls. 
%     tWls               Vector specifying which peak wavelengths are included
%                        in tSpds. Its length must equal the width of tSpds.
%
% Output:
%     predictedProportions   4-column matrix, where each row contains the
%                            predicted proportion for a given row of 
%                            stimParams. Each column represents a different
%                            outcome, in the following order:
%                            1: primary too red/test too bright
%                            2: primary too red/test too dim
%                            3: primary too green/test too bright
%                            4: primary too green/test too dim
%
% Optional key/value pairs
%     'refSpd'         -Spd used as a reference for computing opponent
%                       contrasts. If empty, the opponent contrast of the
%                       primary mixture is calculated relative to the test. 
%                       Default is [].

% History
%     10/21/20  dce  -Adapted from qpPFRM
%     10/23/20  dce  -Edited to not allow complete certainty
%`    11/14/20  dce  -Added option to use without reference spd

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
p.addParameter('refSpd',[],@isnumeric);
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

%% Generate the reference spectrum for opponent contrast calculations
if ~isempty(p.Results.refSpd)
    refLMS = observer.T_cones * p.Results.refSpd;
end 

%% Loop over stimuli and get the proportions for each 
nStim = size(stimParamsVec,1);         % Number of stimulus combinations
predictedProportions = zeros(nStim,4); % Array for storing proportions
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
    
    % Compute opponent contrast of the two lights. If a reference primary
    % is provided, compute each opponent contrast relative to the
    % reference, and return the difference as [test - primary mixture]. 
    % Otherwise, return the contrast of the test relative to the primary mixture.
    if isempty(p.Results.refSpd)
        opponentContrastDiff = LMSToOpponentContrast(observer.colorDiffParams,...
            pLMS,tLMS);
    else
        pOpponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
            refLMS,pLMS);
        tOpponentContrast = LMSToOpponentContrast(observer.colorDiffParams,...
            refLMS,tLMS);
        opponentContrastDiff = tOpponentContrast - pOpponentContrast;
    end
    
    % Compute proportions for luminance judgements: test is not brighter,
    % test is brighter
    pTNbright = normcdf(0,opponentContrastDiff(1),noiseSD);
    pTBright = 1 - pTNbright;  
    
    % Compute proportions for RG judgements: primary is redder, primary is
    % greener
    pPRed = normcdf(0,opponentContrastDiff(2),noiseSD);
    pPGreen = 1 - pPRed;
    
    % Compute proportions for the combined outcomes: too red and too bright, 
    % too red and not too bright, too green and too bright, too green and 
    % not too bright. 
    predictedProportions(ii,:) = [pPRed*pTBright, pPRed*pTNbright,...
        pPGreen*pTBright,pPGreen*pTNbright];
    
    % Don't allow complete certainty
    predictedProportions(ii,predictedProportions(ii,:) > 0.9999) = 0.9999;
    predictedProportions(ii,predictedProportions(ii,:) < 0.0001) = 0.0001;
    predictedProportions(ii,:) = predictedProportions(ii,:)./sum(predictedProportions(ii,:));
    if any(predictedProportions(ii,:)==0 | predictedProportions(ii,:)== 1)
        error('Probability normalization algorithm does not work');
    end 
end