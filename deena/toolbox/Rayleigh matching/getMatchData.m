function [testSpds,primarySpds,testIntensities,primaryRatios] =...
    getMatchData(fName,varargin)
% Helper function that takes in a data file produced by OLRayleighMatch and
% returns the predicted primary and test spds associated with its matches.
% Returns all of the matches by default, but you can also specify to return
% the average spds of all available matches. Also returns the primary ratio
% and test intensity, relative to the scaled OneLight spectra.
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

% Parse input
p = inputParser;
p.addParameter('averageSpds',false,@(x)(islogical(x)));
p.parse(varargin{:});

% Data arrays
% Assume an equal number of matches for each interleaved staircase
trialData = load(fName);   % OLRayleighMatch dataset
spdLength = size(trialData.primarySpds,1);
nMatches = length(trialData.dataArr{1}.matchPositions(:,1));
nInterleaved = length(trialData.dataArr);

testSpds = zeros(spdLength,nMatches,nInterleaved);    % Output test spds
primarySpds = zeros(spdLength,nMatches,nInterleaved); % Output primary spds
primaryRatios = zeros(nMatches,nInterleaved);         % Primary ratio settings
testIntensities = zeros(nMatches,nInterleaved);       % Reference intensity settings

% Find match position indices
for kk = 1:nInterleaved
    if ~isempty(trialData.dataArr{kk}.matchPositions)
        tMatchInds = trialData.dataArr{kk}.matchPositions(:,1);
        pMatchInds = trialData.dataArr{kk}.matchPositions(:,2);
        
        % Define the spd arrays we're searching in
        testArr = trialData.testSpds;
        primaryArr = trialData.primarySpds;
        
        for i = 1:length(tMatchInds)
            % Find spds corresponding to the match position indices. If an index is
            % not an integer, the spd is the average of the spds of the indices
            % above and below it.
            test = mean([testArr(:,ceil(tMatchInds(i))), ...
                testArr(:,floor(tMatchInds(i)))],2);
            testSpds(:,i,kk) = test;
            
            primary = mean([primaryArr(:,ceil(pMatchInds(i))),...
                primaryArr(:,floor(pMatchInds(i)))],2);
            primarySpds(:,i,kk) = primary;
            
            % Find primary ratio and test intensity for each spd. If an index is
            % not an integer, average the scale factors of the indices above and
            % below it.
            primaryRatio = mean([trialData.p1Scales(ceil(pMatchInds(i))),...
                trialData.p1Scales(floor(pMatchInds(i)))]);
            primaryRatios(i,kk) = primaryRatio;
            
            testIntensity = mean([trialData.testScales(ceil(tMatchInds(i))),...
                trialData.testScales(floor(tMatchInds(i)))]);
            testIntensities(i,kk) = testIntensity;
        end
    end
end
% Average if we're looking for averages
if p.Results.averageSpds
    testSpds = mean(testSpds,[2 3]);
    primarySpds = mean(primarySpds,[2 3]);
end
end
