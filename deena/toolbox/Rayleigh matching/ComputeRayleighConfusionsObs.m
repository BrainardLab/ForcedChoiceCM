% Compute locus of confusions in intensity-ratio plot
%
% Syntax:
%    [testIntensity,mixingRatio,matchDiff] =
%    ComputeRayleighConfusions()
%
% Description:
%    Take lambdaMax values and generate receptor fundamentals. Then loop
%    over all test intensities and mixing ratio combinations and compute a
%    measure of color difference between test and corresponding match.
%
%    Many key parameters are specified within this routine rather than
%    passed, because this is a tutorial script.  These include primary
%    wavelengths, matching primary intensities, parameters describing color
%    difference calculation, etc.
%
%    The color difference is computed based on vector length in an
%    post-receptoral contrast space, with different weights applied to the
%    different post-receptoral contrast directions. It is a very rough and
%    ready calculation, but this aspect is not key to demonstrate the
%    principles we are interested in here.
%
% Inputs:
%    observer         - Struct with observer parameters. Must contain the
%                       following fields: colorDiffParams, T_cones. See
%                       genRayleighObserver for details.
%    p1Spd            - Column vector representation of the first primary
%    p2Spd            - Column vector representation of the second primary
%    testSpd          - Column vector representation of the test spd 
%    testIntensityRange       -Row vector of test intensities.  Arbitrary
%                              units.  Values between 0 and 1 are about
%                              right given the way the other parameters are
%                              set.
%    primaryRatioRange        -Row vector of g/r mixing ratios. 0 means all
%                              green primary, 1 means all red. Here green
%                              and red are really defined by the
%                              wavelengths of the two matching primaries
%                              defined in the parameters for this routine.
%
% Outputs:
%    testIntensity             Matrix where entry i,j is the test intensity
%                              given by the ith intensity in testIntensityRange,
%                              and j indexes the mixing ratios.
%    mixingRatio               Matrix where entry i,j is the mixingRatio
%                              given by the jth tentry of mixingRatioRange,
%                              and i indexes the test intensities
%    matchDiff                 Matrix of color differences, where entry i,j
%                              corresponds to the test intensity and mixing
%                              ratio in entry i,j of matrices testIntensity
%                              and mixingRatio.

% History:
%   11/16/20  dce  Adapted from ComputeRayleighConfusions

function [testIntensity,mixingRatio,matchDiff] = ...
    ComputeRayleighConfusionsObs(observer,p1Spd,p2Spd,testSpd,...
    testIntensityRange,mixingRatioRange) 

% I fussed with these to rotate the D line to be horizontal in the plot.
% In real life, they are parameters of the apparatus.
matchLowRawIntensity = 0.1;
overallFactor = 1.65;
relativeHighMatchFactor = 20;
matchIntensity1 = overallFactor*matchLowRawIntensity;
matchIntensity2 = overallFactor*relativeHighMatchFactor*matchIntensity1;

p1Spd = p1Spd.*matchIntensity1;
p2Spd = p2Spd.*matchIntensity2;

% Act like we have an added background that suppresses S cone signals.
% Otherwise small S cone differences explode when we compute contrast,
% because of small denominator.
addedBackgroundCones = [0 0 1]';

% Compute cone respones to test and match
%
% We just do this for all possible test intensities and match ratios, as
% specified in the parameters section.

% Construct each test and compute cone responses
testSpectrum = cell(1,length(testIntensityRange));
testCones = cell(1,length(testIntensityRange));
for ii = 1:length(testIntensityRange)
    testIntensity = testIntensityRange(ii);
    testSpectrum{ii} = testSpd.* testIntensity;
    testCones{ii} = observer.T_cones*(testSpectrum{ii}) + addedBackgroundCones;
end

% Construct each primary mixture and compute cone responses
matchSpectrum = cell(1,length(mixingRatioRange));
matchCones = cell(1,length(mixingRatioRange));
for jj = 1:length(mixingRatioRange)
    mixingRatio = mixingRatioRange(jj);
    matchSpectrum{jj} = (1-mixingRatio)*p1Spd + mixingRatio*p2Spd;
    matchCones{jj} = observer.T_cones*(matchSpectrum{jj})+addedBackgroundCones;
end

% Compute a measure of color difference for each test/match pairing. Uses
% the mixture as the adapting background, and computes the opponent contrast
% of the test relative to it. 
matchDiff = zeros(length(testIntensityRange),length(mixingRatioRange));
mixingRatio = zeros(length(testIntensityRange),length(mixingRatioRange));
testIntensity = zeros(length(testIntensityRange),length(mixingRatioRange));
for ii = 2:length(testIntensityRange)
    for jj = 1:length(mixingRatioRange)
        opponentContrasts = LMSToOpponentContrast(observer.colorDiffParams,...
            testCones{ii},matchCones{jj});
        testIntensity(ii,jj) = testIntensityRange(ii);
        mixingRatio(ii,jj) = mixingRatioRange(jj);
        matchDiff(ii,jj) = norm(opponentContrasts);
    end
end
matchDiff = matchDiff(2:end,:);
mixingRatio = mixingRatio(2:end,:);
testIntensity = testIntensity(2:end,:);
end