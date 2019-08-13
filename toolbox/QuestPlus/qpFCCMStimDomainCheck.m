function paramsOK = qpFCCMStimDomainCheck(stimParamsVec,type,stimParamsStruct)
% Stimulus parameter check for force-choice color matching
%
% Usage:
%     paramsOK = qpFCCMStimDomainCheck(stimParamsVec,type,stimParamsStruct)
%
% Description:
%     Check whether passed parameters are valid stimuli for forcec choice
%     color matching
%
% Inputs:
%     stimParamsVec     Stimulus parameter vec for the FCCM experiment
%     type              String, parameter type
%
% Output:
%     paramsOK          Boolean, true if parameters are OK and false otherwise.

% 07/22/17  dhb  Wrote it.

%% Assume ok
paramsOK = true;

switch (type)
    case 'basic'
        [testParams,comparison1Opponent,comparison2Opponent] = StimulusVecToParams(type,stimParamsVec,stimParamsStruct);
        
        % Only go in one opponent direction in each comparison
        index1 = find(comparison1Opponent == 0);
        if (length(index1) < 2)
            paramsOK = false; return;
        end
        index2 = find(comparison2Opponent == 0);
        if (length(index2) < 2)
            paramsOK = false; return;
        end
        
        index1 = find(comparison1Opponent ~= 0);
        index2 = find(comparison2Opponent ~= 0);
        if (~isempty(index1))
            if (~isempty(index2))
                if (index1 ~= index2)
                    paramsOK = false; return;
                end
            end
        end
        
        % Two comparisons should be in same direction
        if (comparison1Opponent(index1) ~= 0 & comparison1Opponent(index1) ~= 0)
            if (sign(comparison1Opponent(index1)) ~= sign(comparison2Opponent(index2)))
                paramsOK = false; return;
            end
        end
        
        % Don't run trials when two comparisons are same stimulus
        %
        % This checks for the case where all entries are zero and thus
        % the same.
        if (isempty(index1) & isempty(index2))
            paramsOK = false; return;
        end
        
        % If both are non-zero, then check if they are the same.
        if (~isempty(index1) & ~isempty(index2))
            if (comparison1Opponent(index1) == comparison2Opponent(index2))
                paramsOK = false; return;
            end
        end
        
    otherwise
        error('Unknown parameter vector type requested');
        
end


