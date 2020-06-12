function [testSpd, primarySpd] = getSingleMatchData(fName)
% Helper function that takes in a data file produced by OLRayleighMatch and
% returns the predicted primary and test spds associated with the first 
% match. 
%
% Syntax:
%   getSingleMatchData(fName)
%
% Inputs:
%    fName     - Full filename of an OLRayleighMatch data file (.mat)
%
% Outputs:
%    test      - Predicted test spd for the match 
%    primary   - Predicted primary spd for the match 
%
% Optional key-value pairs:
%    None 

% History:
%   06/12/20  dce       Wrote it.

theData = load(fName); 
testSpd = theData.testSpdsPredicted(:,theData.matchPositions(1,1)); 
primarySpd = theData.primarySpdsPredicted(:,theData.matchPositions(1,2));
end 

