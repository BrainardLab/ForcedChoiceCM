function x = StimulusParamsToVec(type,stimulusParams,comparison1Opponent,comparison2Opponent)
% Convert stimulus parameters to vector form
%
% Synopsis:
%   x = StimulusParamsToVec(type,stimulusParams,comparison1Primary,comparison2Primary)
%
% Description:
%   Our goal is to used forced choice color similarity judgments to
%   determine observer parameters.  Sometimes we want those parameters
%   as a vector, and sometimes in more readable form.  This routine goes from
%   the readable form to the vector.
%
% Inputs:
%   type                    - Type of vector to set up.
%                             'basic': 
%   stimulusParams          - Stimulus parameter structure.
%
% Outputs:
%   x                       - Parameters as vector. 
%
% Optional key value pairs:
%   None.
%
% See also:
%   StimulusVecToParams, ComputeCIEConeFundamentals
%

% History:
%   08/09/19  dhb  Wrote it, because I have to do one fun thing this summer.

switch (type)
    case 'basic'
        x = zeros(1,9);
        x(1) = stimulusParams.testParams.testWavelength;
        x(2) =  stimulusParams.matchApparatusParams.primaryWavelength1;
        x(3) =  stimulusParams.matchApparatusParams.primaryWavelength2;
        x(4:6) = comparison1Opponent';
        x(7:9) = comparison2Opponent';
        
    otherwise
        error('Unknown parameter vector type requested');
        
end
