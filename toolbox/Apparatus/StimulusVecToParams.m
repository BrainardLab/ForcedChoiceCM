function [testParams,matchApparatusParams,comparison1Opponent,comparison2Opponent] = StimulusVecToParams(type,x,stimulusParams)
% Convert vector of stimulus parameters to an understandable representation
%
% Synopsis:
%   [testParams,matchApparatusParams,comparison1Opponent,comparison2Opponent] = StimulusVecToParams(type,x,stimulusParams)
%
% Description:
%   Our goal is to used forced choice color similarity judgments to
%   determine observer parameters.  Sometimes we want those parameters
%   as a vector, and sometimes in a more human readable format.  This
%   routine produces the human readable format.
%
% Inputs:
%   type                    - Type of vector to set up.
%                             'basic':
%   x                       - Parameters as vector.
%   params                  - Base stimulus parameter structure.
%
% Outputs:
%   testParams              - Test parameter structure
%   matchApparatusParams    - Apparatus parameters structure
%   comparison1Opponent     - First comparison opponent contrast
%   comparison2Opponent     - Second comparison opponent contrast
%
% Optional key value pairs:
%   None.
%
% See also:
%   StimulusParamsToVec
%

% History:
%   08/09/19  dhb  Wrote it, because I have to do one fun thing this summer.
%   01/02/20  dhb  Added match parameters to this.

% Examples:
%{

%}

switch (type)
    case 'basic'
        stimulusParams.testParams.testWavelength = x(1);
        testParams = SetTestParams(stimulusParams.testParams);
        
        stimulusParams.matchApparatusParams.primaryWavelength1 = x(2);
        stimulusParams.matchApparatusParams.primaryWavelength2 = x(3);
        matchApparatusParams = SetMatchApparatusParams(stimulusParams.matchApparatusParams);
        
        comparison1Opponent = x(4:6)';
        comparison2Opponent = x(7:9)';
        
    otherwise
        error('Unknown parameter vector type requested');
        
end


