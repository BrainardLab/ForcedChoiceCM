% Illustrates how Rayleigh matches vary with changes in cone lambda max and
% optical density
%
% Description:
%   Simulate Rayleigh match performance for variations in L and M cone
%   peak sensitivity and optical density, as described by Thomas and Mollon
%   (2003). Plots Rayleigh match performance in the form of a Pitt diagram.
%
%   The simulated anomaloscope allows adjustment of a monochromatic test
%   and the ratio of two monochromatic primaries in the match. The routine
%   computes the cone responses to the test and match and from these a
%   color difference.  Matches are predicted for test intensity and mixing
%   ratio parameters where the color difference is below a criterion.
%
%   The locus of matches is plotted in a Pitt diagram, where the x-axis is
%   the mixing ratio and the y-axis is the test intensity. The output
%   diagram reproduces the qualitative features of the one that came in the
%   manual for our anamoloscope.
%
%   This routine produces two plots that come close to the Thomas and
%   Mollon results. The first plots Rayleigh match performance for
%   hypothetical dichromats with varying L/M cone lambda maxes, while the
%   second plots Rayleigh matches for dichromats with varying optical
%   densities.
%
%   The color difference model is very simple and is briefly described in
%   the header comments for routine ComputeConfusions, which is at the
%   bottom of this file.

% History
%   07/03/19  dhb       Wrote it.
%   09/03/19  dhb, dce  Modified to use Asano et al. individual difference
%                       parameters, but in the end does the same thing.
%                       However, and enterprising person can now examine
%                       the effect of changing photopigment density.
%   09/12/19  dce       Modified to include specific plots of Thomas and
%                       Mollon results
%% Clear
clear; close all;

%% The first case plots effects of variation in cone lambda max for a 
%% hypothetical dichromat. We begin by defining our parameters

% Cone lambda maxes, taken from Fig. 1 of the Thomas and Mollon paper.
lambdaMaxes1 = [  ...
    [531 531 420.7]' ...
    [541 541 420.7]' ...
    [551 551 420.7]' ...
    [561 561 420.7]'];

% We specify the cones as a shift relative to a nomogram generated lambda
%  max. These base values are set equal to the test values
baseLambdaMaxes1 = [  ...
    [531 531 420.7]' ...
    [541 541 420.7]' ...
    [551 551 420.7]' ...
    [561 561 420.7]'];

% You can also allow the specified photopigment density to
% vary.  Enter these as percent changes relative to nominal
% values. Can be positive or negative. Here, photopigment densities are set
% to 0
dphotopigments1 = zeros(3,4);

% Plotting parameters
colors1 = [ 'r' 'g' 'b' 'y'];
legend1 = {'531' '541' '551' '561'};
title1 = 'Lambda Max Variation'; 
testIntensityRange1 = 0:0.001:0.35; 
plotMatches(lambdaMaxes1, baseLambdaMaxes1, dphotopigments1, colors1, legend1, title1, testIntensityRange1)

%% The second case plots effects of optical density variation for a  
%% hypothetical dichromat

% All lambda max and base lambda max values are the same
lambdaMaxes2 = [[561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]'];

baseLambdaMaxes2 = [[561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]' ...
    [561 561 420.7]'];

% Optical density shifts are taken from Fig. 1 of the Thomas and Mollon
% paper: 0.05, 0.2, 0.35, 0.5, 0.65. The values here are calculated as 
% percent changes from the average optical density of 0.3.
dphotopigments2 = [  ...
    [-83.33 -83.33 0]' ...
    [-33.33 -33.33 0]' ...
    [16.67 16.67 0]' ...
    [66.67 66.67 0]'...
    [116.67 116.67 0]'];

colors2 = [ 'r' 'g' 'b' 'y' 'm' ];
legend2 = { '-83.33%' '-33.33%' '16.67%' '66.67%' '116.67%' };
title2 = 'Optical Density Variation'; 
testIntensityRange2 = 0.:0.001:0.35; 
plotMatches(lambdaMaxes2, baseLambdaMaxes2, dphotopigments2, colors2, legend2, title2, testIntensityRange2)

function plotMatches(lambdaMaxes, baseLambdaMaxes, dphotopigments, colors, theLegend, theTitle, testIntensityRange)
% Convert specified lambda max values to shifts from the nominal CIE
% standard values.
for ii = 1:size(lambdaMaxes,2)
    indDiffParams(ii).dlens= 0;
    indDiffParams(ii).dmac = 0;
    indDiffParams(ii).dphotopigment = dphotopigments(:,ii)';
    indDiffParams(ii).lambdaMaxShift = lambdaMaxes(:,ii)' - baseLambdaMaxes(:,ii)';
    indDiffParams(ii).shiftType = 'linear';
end

% Threshold difference below which it is a match
% Fussed with this by hand to adjust plot to taste.
thresholdVal = 0.12;

% Apparatus range parameters.
% Mixing ratio of zero is all green primary, 1 is all red
% (assuming shorter primary is first one specified in routine
% below.)
mixingRatioRange = 0:0.001:1;

%% Loop to calculate matching locus for each set of cones
%
% For each set of specified cone pigment lambda max, this
% computes the matching range and adds to the Pitt diagram
% plot.
theFigure = figure; clf; hold on
for kk = 1:size(lambdaMaxes,2)
    % Function below does the work, based on lambdaMax.
    % Most operating parameters are set in the function itself.
    [testIntensity{kk},mixingRatio{kk},matchDiff{kk}] = ComputeConfusions(baseLambdaMaxes(:,kk),indDiffParams(kk),testIntensityRange,mixingRatioRange);
    
    % This plot will show the color difference as a function of mixing ratio
    % and test intensity, one plot per set of lambda max values.  I found these
    % useful for development but not all that instructive in the end, so
    % they conditional and off by default.
    diffPlots = false;
    if (diffPlots)
        figure; clf; hold on
        mesh(mixingRatio{kk},testIntensity{kk},matchDiff{kk});
        colormap(winter)
        view([2 14]);
        zlim([0 2]);
        xlim([min(mixingRatioRange) max(mixingRatioRange)]);
        ylim([min(testIntensityRange) max(testIntensityRange)]);
        xlabel(' Mixing Ratio (0 -> green; 1 -> red)');
        ylabel('Test Intensity');
        zlabel('Color Difference');
    end
    
    figure(theFigure);
    index = find(matchDiff{kk} < thresholdVal);
    plot(mixingRatio{kk}(index),testIntensity{kk}(index),[colors(kk) 'o'],'MarkerFaceColor',colors(kk));
    
end

% Finish off the plot
figure(theFigure);
xlim([min(mixingRatioRange) max(mixingRatioRange)]);
ylim([min(testIntensityRange) max(testIntensityRange)]);
xlabel(' Mixing Ratio (0 -> green; 1 -> red)');
ylabel('Test Intensity');
axis('square')
legend(theLegend);
title(theTitle)
end

% Compute locus of confusions in intensity-ratio plot
%
% Syntax:
%    [testIntensity,mixingRatio,matchDiff] = ComputeConfusions(lambdaMax,indDiffParams,testIntensityRange,mixingRatioRange)
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

function [testIntensity,mixingRatio,matchDiff] = ComputeConfusions(lambdaMax,indDiffParams,testIntensityRange,mixingRatioRange)
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
