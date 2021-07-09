function [testSpds,primarySpds,testIntensities,primaryRatios] =...
    getMatchData_old(fName,varargin)
% This is an old version of the function getMatchData, for use with older
% OLRayleighMatch files which were saved in a different format. Called by 
% OLAnalyzeRayleighMatch_old. See getMatchData for full description. 
%
% Syntax:
%   getMatchData(fName)
%
% Inputs:
%    fName     - Full filename of an OLRayleighMatch data file (.mat)
%
% Outputs:
%    testSpds        - Predicted test spd(s) for the match
%    primarySpds     - Predicted primary spd(s) for the match
%    testIntensities - Test intensity settings at the match
%    primaryRatios   - Primary ratio settings at the match
%
% Optional key-value pairs:
%    averageSpds - Logical indicating to return the average of all spds in
%                  the file (note that test intensities and primary ratios
%                  are not averaged). Default is false.

% History:
%   06/12/20  dce       Wrote it.
%   06/16/20  dce       Modified to return one or multiple spd pairs, fixed
%                       for non-integer match positions
%   07/06/20  dce       Added primary ratio and test intensity outputs.
%   07/17/20  dce       Added option to use nominal spds, not predicted
%   07/24/20  dce       Added handling of case where no matches were made
%   08/07/20  dce       Added spd averaging option
%   06/02/21  dce       Edited to reflect changes to OLRayleighMatch file
%                       structure
%   06/22/21  dce       Edited to calculate matches based on last
%                       reversals, not last settings

% Parse input
p = inputParser;
p.addParameter('averageSpds',false,@(x)(islogical(x)));
p.parse(varargin{:});

% Data arrays
% Assume an equal number of matches for each interleaved staircase
trialData = load(fName);   % OLRayleighMatch dataset
spdLength = size(trialData.primarySpdsPredicted,1);
nMatches = length(trialData.matchPositions(:,1));
nInterleaved = 1;

testSpds = zeros(spdLength,nMatches,nInterleaved);    % Output test spds
primarySpds = zeros(spdLength,nMatches,nInterleaved); % Output primary spds
primaryRatios = zeros(nMatches,nInterleaved);         % Primary ratio settings
testIntensities = zeros(nMatches,nInterleaved);       % Reference intensity settings

% Find match position indices
for kk = 1:nInterleaved
     if ~isempty(trialData.matchPositions)
        tMatchInds = trialData.matchPositions(:,1);
        pMatchInds = trialData.matchPositions(:,2);
        
        % Define the spd arrays we're searching in
        testArr = trialData.testSpdsPredicted;
        primaryArr = trialData.primarySpdsPredicted;
        
        for i = 1:length(tMatchInds)
            % Find spds corresponding to the match position indices. If an index is
            % not an integer, the spd is a weighted average of the spds of the indices
            % above and below it.
            testScale = (1-(ceil(tMatchInds(i))-tMatchInds(i)));
            test = ...
                (testArr(:,ceil(tMatchInds(i)))*testScale+...
                testArr(:,floor(tMatchInds(i)))*(1-testScale));
            testSpds(:,i,kk) = test;
            
            pScale = (1-(ceil(pMatchInds(i))-pMatchInds(i)));
            primary = ...
                (primaryArr(:,ceil(pMatchInds(i)))*pScale+...
                primaryArr(:,floor(pMatchInds(i)))*(1-pScale));
            primarySpds(:,i,kk) = primary;
            
            % Find primary ratio and test intensity for each spd. If an index is
            % not an integer, take the weighted average of the indices
            % above and below it.
            primaryRatios(i,kk) = trialData.matches(i,2);            
            testIntensities(i,kk) = trialData.matches(i,1);
        end
    end
end
% Average spds if desired 
if p.Results.averageSpds
    testSpds = mean(testSpds,[2 3]);
    primarySpds = mean(primarySpds,[2 3]);
end
end