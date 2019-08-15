function [predictedProportions] = qpPFFCCM(stimParamsVec,psiParamsVec,S,stimVectorType,stimParamsStruct,psiVectorType,psiParamsStruct,psiParamsStructRef,adaptationSpd,varargin)
% Psychometric function for forced choice color matching
%
% Usage:
%     [predictedProportions] = qpPFFCCM(stimParamsVec,psiParamsVec,S,stimVectorType,stimParamsStruct,psiVectorType,psiParamsStruct,psiParamsStructRef,adaptationSpd)
%
% Description:
%     Psychometric function for forced choice color matching.  
%
% Input:
%     stimParamsVec          Matrix, with each row being a vector of stimulus parameters.
%     psiParamsVec           Row vector of parameters
%     S                      Wavelength sampling
%     stimVectorType         String specifying type of stimulus parameter vector
%     stimParamsStruct       Full stimulus parameters struct
%     psiVectorType          String specifying vector type of psi parameter vector
%     psiParamsStruct        Full observer parameters structure
%     psiParamsStructRef     Full observer parameters structure for a reference observer.
%     adaptationSpd          Spectral power distribution of adapting light
%
% Output:
%     predictedProportions   Matrix, where each row is a vector of predicted proportions
%                            for each outcome.
%                              First entry of each row is for first comparison chosen (outcome == 1)
%                              Second entry of each row is for second comparison chosen (outcome == 2)
%
% Optional key/value pairs
%     None.

% 07/03/17  dhb  Wrote it

%% Parse input
p = inputParser;
p.addRequired('stimParamsVec',@isnumeric);
p.addRequired('psiParamsVec',@isnumeric);
p.addRequired('S',@isnumeric);
p.addRequired('stimVectorType',@ischar);
p.addRequired('stimParamsStruct',@isstruct);
p.addRequired('psiVectorType',@ischar);
p.addRequired('psiParamsStruct',@isstruct);
p.addRequired('psiParamsStructRef',@isstruct);
p.addRequired('adaptationSpd',@isnumeric);
%p.addParameter('slope',3,@isscalar);
p.parse(stimParamsVec,psiParamsVec,S,stimVectorType,stimParamsStruct,psiVectorType,psiParamsStruct,psiParamsStructRef,adaptationSpd,varargin{:});

%% Get observer parameters structure from psiParams vector
psiParamsStruct = ObserverVecToParams(psiVectorType,psiParamsVec,psiParamsStruct);

%% Get cone funadmentals
TRef = ComputeObserverFundamentals(psiParamsStructRef.coneParams,S);
T = ComputeObserverFundamentals(psiParamsStruct.coneParams,S);

%% LMS coordinates of adapting light
%
% Here we just care about these with respect to the 
% reference observer, so that we can keep our stimulus
% transformations independent of the current observer 
% parameters.
adaptationLMSRef = TRef*adaptationSpd;
adaptationLMS = T*adaptationSpd;

% Matrix for opponent transformation
% Be careful not to search over opponent space
M = GetOpponentContrastMatrix(psiParamsStructRef.colorDiffParams);

%% Loop over stimuli and get the proportions
nStim = size(stimParamsVec,1);
predictedProportions = zeros(nStim,2);
for ii = 1:nStim
    
    % Get stimulus parameters from stimParams vector
    [testParams,comparison1Opponent,comparison2Opponent] = StimulusVecToParams(stimVectorType,stimParamsVec(ii,:),stimParamsStruct);
    
    % Get reference spectrum.
    referenceSpd = adaptationSpd + testParams.testIntensity*testParams.unitTestSpectrum;
    referenceLMS = T*referenceSpd;
        
    % Get comparison spectra.  We do this with respect to the reference
    % observer, because what we are doing here is converting stimulus
    % opponent representation, which is convenient to think about, into
    % stimulus spectra. We don't want to do this with respect to a
    % floating observer specification.
    referenceLMSRef = TRef*referenceSpd;
    referenceOpponentRef = LMSToOpponentContrast(M,adaptationLMSRef,referenceLMSRef);
    
    % Find comparison primaries
    comparison1LMSRef = OpponentContrastToLMS(M,adaptationLMSRef,referenceOpponentRef+comparison1Opponent);
    comparison2LMSRef = OpponentContrastToLMS(M,adaptationLMSRef,referenceOpponentRef+comparison2Opponent);
    
    [comparison1Primary,comparison1Spd] = FindMetamer(stimParamsStruct.matchApparatusParams,TRef,comparison1LMSRef);
    [comparison2Primary,comparison2Spd] = FindMetamer(stimParamsStruct.matchApparatusParams,TRef,comparison2LMSRef);
    comparison1LMS = T*comparison1Spd;
    comparison2LMS = T*comparison2Spd;
    
    % Do the main work
    predictedProportions(ii,1) = ...
        ComputeChoiceLikelihood(psiParamsStruct,M,adaptationLMS,referenceLMS,comparison1LMS,comparison2LMS);
    
    % Fill in complement
    predictedProportions(ii,2) = 1-predictedProportions(ii,1);  
end

%% Don't allow complete certainty
predictedProportions(predictedProportions > 0.9999) = 0.9999;
predictedProportions(predictedProportions < 0.0001) = 0.0001;

