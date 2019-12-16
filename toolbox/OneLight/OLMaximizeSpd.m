function [spdMax,primaryScale] = OLMaximizeSpd(cal,spdIncrRel,varargin)


%% Parse input
%
% This currently sets the primary and test wavelengths.
p = inputParser;
p.addParameter('minScaleFactor', 0.1,@(x) (isnumeric(x)));
p.addParameter('maxScaleFactor', 10, @(x) (isnumeric(x)));
p.addParameter('errorTolerance', 0.05, @(x) (isnumeric(x)));
p.addParameter('nScaleFactors', 100, @(x) (isnumeric(x)));
p.addParameter('lambda', 100, @(x) (isnumeric(x)));
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
primaryScales = linspace(p.Results.minScaleFactor,p.Results.maxScaleFactor,p.Results.nScaleFactors);
for ii = 1:length(primaryScales)
    primaryScale = primaryScales(ii);
    targetSpd = primaryScale*spdIncrRel+darkSpd;
    [targetPrimary,targetSpdPredicted] = OLSpdToPrimary(cal, targetSpd, 'lambda', p.Results.lambda, ...
        'whichSpdToPrimaryMin', 'leastSquares', ...
        'verbose', false);
    relativeRMSEError(ii) = sqrt(sum( ((targetSpdPredicted-targetSpd)/mean(targetSpd)).^2 ) );
end

%% Find minimum error to determine a good maximum scale factor
minVal = min(relativeRMSEError);

% Find the largest scale factor where erros is less than tolerance
% more than max
targetRMSEVal = (1 + p.Results.errorTolerance)*minVal;
indices = find(relativeRMSEError < targetRMSEVal);
candidateScales = primaryScales(indices);
primaryScale = max(candidateScales);
spdMax = primaryScale*spdIncrRel;