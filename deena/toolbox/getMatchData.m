function [testSpds,primarySpds,primaryRatio,testIntensity] = getMatchData(fName, varargin)
% Helper function that takes in a data file produced by OLRayleighMatch and
% returns the predicted primary and test spds associated with its matches.
% Returns all of the matches by default, but you can also specify an index
% so it returns a specific match. Also returns the primary ratio and test
% intensity, relative to the base OneLight spectra
%
% Syntax:
%   getMatchData(fName)
%
% Inputs:
%    fName     - Full filename of an OLRayleighMatch data file (.mat)
%
% Outputs:
%    test      - Predicted test spd(s) for the match
%    primary   - Predicted primary spd(s) for the match
%
% Optional key-value pairs:
%    ind       - Positive integer index of the desired match. If an index
%                is specified, the function returns just one pair of spds
%                instead of all available matches. The default value is 0,
%                which is ignored.

% History:
%   06/12/20  dce       Wrote it.
%   06/16/20  dce       Modified to return one or multiple spd pairs, fixed
%                       for non-integer match positions
%   07/06/20  dce       Added primary ratio and test intensity outputs.

p = inputParser;
p.addParameter('ind',0,@(x)(isnumeric(x)));
p.parse(varargin{:});
ind = p.Results.ind;

theData = load(fName);   % OLRayleighMatch dataset
testSpds = [];           % Output test spds
primarySpds = [];        % Output primary spds
primaryRatios = [];      % Primary ratio settings
testIntensity = [];      % Test intensity settings

% Find match position indices
if ind == 0   % Return spds for all matches
    tMatchInds = theData.matchPositions(:,1);
    pMatchInds = theData.matchPositions(:,2);
else          % Return spds for a single match
    tMatchInds = theData.matchPositions(ind,1);
    pMatchInds = theData.matchPositions(ind,2);
end


for i = 1:length(tMatchInds)
    % Find spds corresponding to the match position indices. If an index is 
    % not an integer, the spd is the average of the spds of the indices 
    % above and below it.
    test = mean([theData.testSpdsPredicted(:,ceil(tMatchInds(i))),...
        theData.testSpdsPredicted(:,floor(tMatchInds(i)))],2);
    testSpds = [testSpds, test];
    
    primary = mean([theData.primarySpdsPredicted(:,ceil(pMatchInds(i))),...
        theData.primarySpdsPredicted(:,floor(pMatchInds(i)))],2);
    primarySpds = [primarySpds, primary];
    
    % Find primary ratio and test intensity for each spd. If an index is
    % not an integer, average the scale factors of the indices above and
    % below it.
    p1Ratio = mean(theData.p1Scales(ceil(pMatchInds(i))),...
        theData.p1Scales(floor(pMatchInds(i)))); 
    testIntensity = mean(theData.testScales(ceil(tMatchInds(i))),...
        theData.testScales(floor(tMatchInds(i)))); 
end




end

