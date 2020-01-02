function [predictedProportions] = qpPFFCCM(stimParamsVec,psiParamsVec,S,stimVectorType,stimParamsStruct,psiVectorType,psiParamsStruct,psiParamsStructRef,adaptationSpd,TRef,adaptationLMSRef,varargin)
% Psychometric function for forced choice color matching
%
% Usage:
%     [predictedProportions] = qpPFFCCM(stimParamsVec,psiParamsVec,S,stimVectorType,stimParamsStruct,psiVectorType,psiParamsStruct,psiParamsStructRef,adaptationSpd,TRef,adaptationLMSRef)
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
%     TRef                   Reference color matching functions.  We could
%                            compute these here, but passing speeds things up a little.
%     adaptationLMSRef       adaptation LMS under reference T (TRef).  We
%                            could compute this here, but passing speeds things up a little.
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
% 01/02/20  dhb  Allow match parameters to vary

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
T = ComputeObserverFundamentals(psiParamsStruct.coneParams,S);

%% LMS coordinates of adapting light
%
% Here we just care about these with respect to the 
% reference observer, so that we can keep our stimulus
% transformations independent of the current observer 
% parameters.
adaptationLMS = T*adaptationSpd;

%% Loop over stimuli and get the proportions
nStim = size(stimParamsVec,1);
predictedProportions = zeros(nStim,2);
for ii = 1:nStim
    
    % Get stimulus parameters from stimParams vector
    [testParams,matchApparatusParams,comparison1Opponent,comparison2Opponent] = StimulusVecToParams(stimVectorType,stimParamsVec(ii,:),stimParamsStruct);
    
    % Get reference spectrum.
    referenceSpd = adaptationSpd + testParams.testIntensity*testParams.unitTestSpectrum;
    referenceLMS = T*referenceSpd;
        
    % Get comparison spectra.  We do this with respect to the reference
    % observer, because what we are doing here is converting stimulus
    % opponent representation, which is convenient to think about, into
    % stimulus spectra. We don't want to do this with respect to a
    % floating observer specification.
    referenceLMSRef = TRef*referenceSpd;
    referenceOpponentRef = LMSToOpponentContrast(psiParamsStructRef.colorDiffParams,adaptationLMSRef,referenceLMSRef);
    
    % Find comparison primaries
    comparison1LMSRef = OpponentContrastToLMS(psiParamsStructRef.colorDiffParams,adaptationLMSRef,referenceOpponentRef+comparison1Opponent);
    comparison2LMSRef = OpponentContrastToLMS(psiParamsStructRef.colorDiffParams,adaptationLMSRef,referenceOpponentRef+comparison2Opponent);
    
    [comparison1Primary,comparison1Spd] = FindMetamer(matchApparatusParams,TRef,comparison1LMSRef);
    [comparison2Primary,comparison2Spd] = FindMetamer(matchApparatusParams,TRef,comparison2LMSRef);
    comparison1LMS = T*comparison1Spd;
    comparison2LMS = T*comparison2Spd;
    
    % Do the main work
    predictedProportions(ii,1) = ...
        ComputeChoiceLikelihood(psiParamsStruct,adaptationLMS,referenceLMS,comparison1LMS,comparison2LMS);
    
    % Fill in complement
    predictedProportions(ii,2) = 1-predictedProportions(ii,1);  
end

%% Don't allow complete certainty
predictedProportions(predictedProportions > 0.9999) = 0.9999;
predictedProportions(predictedProportions < 0.0001) = 0.0001;

