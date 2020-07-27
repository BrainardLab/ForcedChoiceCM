function [testSpds,primarySpds,testIntensities,primaryRatios] =...
    getMatchData(fName,varargin)
% Helper function that takes in a data file produced by OLRayleighMatch and
% returns the predicted primary and test spds associated with its matches.
% Returns all of the matches by default, but you can also specify an index
% so it returns a specific match. Also returns the primary ratio and test
% intensity, relative to the scaled OneLight spectra.
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
%    nominal   - Logical indicating to return spds from the nominal spds
%                array, not the predicted array. Default is false

% History:
%   06/12/20  dce       Wrote it.
%   06/16/20  dce       Modified to return one or multiple spd pairs, fixed
%                       for non-integer match positions
%   07/06/20  dce       Added primary ratio and test intensity outputs.
%   07/17/20  dce       Added option to use nominal spds, not predicted
%   07/24/20  dce       Added handling of case where no matches were made

% Example

% Parse input
p = inputParser;
p.addParameter('ind',0,@(x)(isnumeric(x)));
p.addParameter('nominal',false,@(x)(islogical(x)));
p.parse(varargin{:});
ind = p.Results.ind;

% Data arrays
theData = load(fName);   % OLRayleighMatch dataset
testSpds = [];           % Output test spds
primarySpds = [];        % Output primary spds
primaryRatios = [];      % Primary ratio settings
testIntensities = [];      % Test intensity settings

% Find match position indices
if ~isempty(theData.matchPositions)
    if ind == 0   % Return spds for all matches
        tMatchInds = theData.matchPositions(:,1);
        pMatchInds = theData.matchPositions(:,2);
    else          % Return spds for a single match
        tMatchInds = theData.matchPositions(ind,1);
        pMatchInds = theData.matchPositions(ind,2);
    end
    
    % Define the spd arrays we're searching in
    if p.Results.nominal
        testArr = theData.testSpdsNominal;
        primaryArr = theData.primarySpdsNominal;
    else
        testArr = theData.testSpdsPredicted;
        primaryArr = theData.primarySpdsPredicted;
    end
    
    for i = 1:length(tMatchInds)
        % Find spds corresponding to the match position indices. If an index is
        % not an integer, the spd is the average of the spds of the indices
        % above and below it.
        test = mean([testArr(:,ceil(tMatchInds(i))), ...
            testArr(:,floor(tMatchInds(i)))],2);
        testSpds = [testSpds,test];
        
        primary = mean([primaryArr(:,ceil(pMatchInds(i))),...
            primaryArr(:,floor(pMatchInds(i)))],2);
        primarySpds = [primarySpds,primary];
        
        % Find primary ratio and test intensity for each spd. If an index is
        % not an integer, average the scale factors of the indices above and
        % below it.
        primaryRatio = mean([theData.p1Scales(ceil(pMatchInds(i))),...
            theData.p1Scales(floor(pMatchInds(i)))]);
        primaryRatios = [primaryRatios,primaryRatio];
        
        testIntensity = mean([theData.testScales(ceil(tMatchInds(i))),...
            theData.testScales(floor(tMatchInds(i)))]);
        testIntensities = [testIntensities,testIntensity];
    end
end
end

