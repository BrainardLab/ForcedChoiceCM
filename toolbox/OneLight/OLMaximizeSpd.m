function [maxSpdIncrRel,spdScaleFactor] = OLMaximizeSpd(cal,spdIncrRel,varargin)
% Maximize an incremental spectral power distribution, within OneLight gamut
%
% Synopsis:
%    [maxSpdIncrRel,spdScaleFactor] = OLMaximizeSpd(cal,spdIncrRel)
%
% Description:
%     Take in an inremental spd (to be added to OL dark spd) and figure out
%     how it needs to be scaled to be at the upper edge of the gamut.
%
%     It does this by searching over a number of possible scale factors,
%     computed linearly over the range provided (see key/value pairs).
%
%     The routine will allow a bit of maxed out primary values, as long as a
%     measure of fit of predicted versus target spectrum does not get too large.  You
%     can control how large that is relative to the error of a fully in-gamut
%     choice (again, see key/value pairs).
%
%     This routine does assume that there is some scale factor within the searched
%     range that will bring the requested primary into gamut, and does not check whether
%     negative primary values would be required to do so.
%
% Inputs:
%     cal             - One Light calibration structure
%     spdIncrRel      _ The desired incremental spectrum, in the same units
%                       and wavelength spacing as the calibration data.
%
% Outputs:
%     maxSpdIncrRel   - The determined maximum incremental spd, with same
%                       spectral shape as the input.  Add the One Light
%                       dark light spd to this to get the target spd you should
%                       ask the device to produce.
%     spdScaleFactor  - The scale factor applied to the input spd to get
%                       the output spd.
%     darkSpd         - The OneLight dark light spd, extracted from the
%                       calibration structure and returned for convenience.
% 
% Optional key/value pairs:
%   'minScaleFactor'  - Minimum scale factor to apply to passed relative incremental spd.
%                       Default 0.01.
%   'maxScaleFactor'  - Maximum scale factor to apply to passed relative incremental spd.
%                       Default 10.
%   'nScaleFactors'   - How many scale factors to be linspaced between min and max.
%                       Default 200.
%   'errorTolerance'  - Fractional error increase to tolerate after primaries hit 1
%                       Default 0.05.
%   'nKeepComputing'  - How many scale factors to compute for after primaries hit 1.
%                       Default 4.
%   'lambda'          - The lambda parameter passed into OLSpdToPrimary.
%                       Default 0.001.
%   'showPlot'        - Make a plot that helps see what this is doing.
%                       Default false.

% History
%   12/18/20  dhb  Got it working reasonably




%% Parse input
%
% This currently sets the primary and test wavelengths.
p = inputParser;
p.addParameter('minScaleFactor', 0.01, @(x) (isnumeric(x)));
p.addParameter('maxScaleFactor', 10, @(x) (isnumeric(x)));
p.addParameter('nScaleFactors', 200, @(x) (isnumeric(x)));
p.addParameter('errorTolerance', 0.05, @(x) (isnumeric(x)));
p.addParameter('nKeepComputing', 4, @(x) (isnumeric(x)));
p.addParameter('lambda', 0.001, @(x) (isnumeric(x)));
p.addParameter('showPlot', false, @(x) (islogical(x)));
p.parse(varargin{:});

%% Get the "dark" spd, that which comes out when primarys are at 0
darkSpd = OLPrimaryToSpd(cal,zeros(size(cal.computed.D,2),1));

%% Find primaries for a number scale factors
%
% Compute relative error for each.  Relative error
% is used because otherwise it will increase with
% scale factor even if the shape is the same.
%
% Differential mode looks at the effect of the primaries
% relative to the light that comes out when the input
% is set to 0. For our purposes, this is what we want.
%
% Preallocate space
spdScaleFactors = linspace(p.Results.minScaleFactor,p.Results.maxScaleFactor,p.Results.nScaleFactors);
relativeRMSEError = zeros(size(spdScaleFactors));
primaryOutOfRanges = zeros(size(spdScaleFactors));
targetSpds = zeros(length(spdIncrRel),length(spdScaleFactors));
targetSpdsPredicted = zeros(size(targetSpds));
targetPrimaries = zeros(size(cal.computed.D,2),length(spdScaleFactors));

% Loop over scale factors
for ii = 1:length(spdScaleFactors)
    spdScaleFactor = spdScaleFactors(ii);
    targetSpds(:,ii) = spdScaleFactor*spdIncrRel+darkSpd;
    [targetPrimaries(:,ii),targetSpdsPredicted(:,ii)] = OLSpdToPrimary(cal, targetSpds(:,ii), 'lambda', p.Results.lambda, ...
        'whichSpdToPrimaryMin', 'leastSquares', ...
        'verbose', false);
    
    % Compute relative error
    relativeRMSEError(ii) = sqrt(sum( ((targetSpdsPredicted(:,ii)-targetSpds(:,ii))/mean(targetSpds(:,ii))).^2 ) );

    % Keep track of primaries out of range
    if (any(targetPrimaries(:,ii) >= 1-1e-6))
        primaryOutOfRanges(ii) = 1;
    else
        primaryOutOfRanges(ii) = 0;
    end
    
    % Don't bother computing too far once we're clearly out of range
    if (sum(primaryOutOfRanges) > p.Results.nKeepComputing)
        primaryOutOfRanges((ii+1):end) = true;
        relativeRMSEError((ii+1):end) = Inf;
        break;
    end
end

%% Find reasonable value for RMSE error to determine a good maximum scale factor
indices = find(~primaryOutOfRanges);
reasonableRMSEVal = max(relativeRMSEError(indices));

%% Find the largest scale factor where erros is less than tolerance
% more than the reasonable value.
targetRMSEVal = (1 + p.Results.errorTolerance)*reasonableRMSEVal;
indices = find(relativeRMSEError < targetRMSEVal);
candidateScales = spdScaleFactors(indices);
spdScaleFactor = max(candidateScales);
maxSpdIncrRel = spdScaleFactor*spdIncrRel;

%% Diagonostic figure
if (p.Results.showPlot)
    figure; clf; hold on
    plot(spdScaleFactors,relativeRMSEError);
    index = find(spdScaleFactor == spdScaleFactors);
    plot(spdScaleFactor,relativeRMSEError(index),'ro','MarkerFaceColor','r');
    xlabel('Scale Factor');
    ylabel('Relative RMSE Value');
end

