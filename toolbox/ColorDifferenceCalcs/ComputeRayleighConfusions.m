% Compute locus of confusions in intensity-ratio plot
%
% Syntax:
%    [testIntensity,mixingRatio,matchDiff] =
%    ComputeRayleighConfusions(lambdaMax,indDiffParams,testIntensityRange,mixingRatioRange)
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
%    lambdaMax                 Column vector of three receptor photopigment lambda
%                              max (wavelength of peak sensitivity) values, in nm.
%    indDiffParams             Passed to ComputeCIEConeFundamentals.
%                              Ignored if empty.  If you pass this
%                              structure, then lambdaMax should be empty,
%                              and vice-versa.  That is, only adjust the
%                              fundamentals using one of the two available
%                              methods.
%    testIntensityRange        Row vector of test intensities.  Arbitrary
%                              units.  Values between 0 and 1 are about
%                              right given the way the other parameters are
%                              set.
%    mixingRatioRange          Row vector of g/r mixing ratios. 0 means all
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
%   07/04/19  dhb  Made this its own routine.
%   10/10/19  dce  Made separate routine

function [testIntensity,mixingRatio,matchDiff] = ComputeRayleighConfusions(lambdaMax,indDiffParams,testIntensityRange,mixingRatioRange)
% Observer parameters
fieldSizeDegs = 2;
observerAge = 32;
pupilDiameterMM = 3;

% Wavelength sampling. Life is easiest at 1 nm sampling.
S = [380 1 401];
wls = SToWls(S);

% Apparatus parameters.  These match the Nagel in wavelengths.
testWavelength = 590;
matchWavelength1 = 545;
matchWavelength2 = 679;

% I fussed with these to rotate the D line to be horizontal in the plot.
% In real life, they are parameters of the apparatus.
matchIntensity1 = 0.1;
matchIntensity2 = 3;

% Compute indices so that we can set spectra below
testIndex = find(wls == testWavelength);
matchIndex1 = find(wls == matchWavelength1);
matchIndex2 = find(wls == matchWavelength2);

% Color difference computation parameters.
% I fussed with these to make the uncertainty
% regions look a bit like those in our device's
% diagram.
LMRatio = 2;
lumWeight = 4;
rgWeight = 2;
sWeight = 0.5;

% Act like we have an added background that suppresses S cone signals.
% Otherwise small S cone differences explode when we compute contrast,
% because of small denominator.
addedBackgroundCones = [0 0 1]';

% Generate match spectra before application of mixing ratio
matchSpectrum1 = zeros(size(wls)); matchSpectrum1(matchIndex1) = matchIntensity1;
matchSpectrum2 = zeros(size(wls)); matchSpectrum2(matchIndex2) = matchIntensity2;

% Generate the cones
%
% The weird looking call around the CompueCIEConeFundamentals has the net
% effect of putting the cone fundamentals into energy units, and then we
% normalize each to a peak of one.
%
% See ComputeCIEConeFundamentals for more info, and for other ways to shift
% individual difference parameters.
T_cones = EnergyToQuanta(S, ...
    ComputeCIEConeFundamentals(S,fieldSizeDegs,observerAge,pupilDiameterMM,lambdaMax, ...
    [],[],[],[],[],indDiffParams)')';

for ii = 1:size(T_cones,1)
    T_cones(ii,:) = T_cones(ii,:)/max(T_cones(ii,:));
end

% Make diagnostic plot of cone fundamentals?
% FUNDAMENTAL_PLOTS = ;
% figure; clf; hold on;
% plot(SToWls(S),T_cones(1,:),'r','LineWidth',2);
% plot(SToWls(S),T_cones(2,:),'g','LineWidth',2);
% plot(SToWls(S),T_cones(3,:),'b','LineWidth',2);
% xlabel('Wavelength');
% ylabel('Fundamental');

% Compute cone respones to test and match
%
% We just do this for all possible test intensities and match ratios, as
% specified in the parameters section.

% Construct each test and compute cone responses
for ii = 1:length(testIntensityRange)
    testIntensity = testIntensityRange(ii);
    testSpectrum{ii} = zeros(size(wls)); testSpectrum{ii}(testIndex) = testIntensity;
    testCones{ii} = T_cones*testSpectrum{ii};
end

% Construct each match and compute cone responses
for jj = 1:length(mixingRatioRange)
    mixingRatio = mixingRatioRange(jj);
    matchSpectrum{jj} = (1-mixingRatio)*matchSpectrum1 + (mixingRatio)*matchSpectrum2;
    matchCones{jj} = T_cones*matchSpectrum{jj};
end

% Compute a measure of color difference for each test/match pairing
%
% We'll take the test as contributing to the adapting background and compute difference as
% cone contrast with respect to that plus the added background as specfied
% above.
for ii = 1:length(testIntensityRange)
    for jj = 1:length(mixingRatioRange)
        effectiveBackgroundCones{ii} = testCones{ii} + addedBackgroundCones;
        coneContrastDiff = (testCones{ii}-matchCones{jj})./effectiveBackgroundCones{ii};
        
        % Approximate three post-receptoral constrasts
        lumContrast(ii,jj) = (LMRatio*coneContrastDiff(1)+coneContrastDiff(2))/(LMRatio+1);
        rgContrast(ii,jj) = coneContrastDiff(1)-coneContrastDiff(2);
        sContrast(ii,jj) = coneContrastDiff(3);
        
        % Take weighted sum of squares.  I'm making weights up on grounds
        % that rg is most sensitive, lum next, and s last.  Very back of
        % envelope and may not be right for uniform fields.
        testIntensity(ii,jj) = testIntensityRange(ii);
        mixingRatio(ii,jj) = mixingRatioRange(jj);
        matchDiff(ii,jj) = sqrt((lumWeight*lumContrast(ii,jj))^2 + (rgWeight*rgContrast(ii,jj))^2 + (sWeight*sContrast(ii,jj))^2);
    end
end

end