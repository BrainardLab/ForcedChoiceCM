function spd = makeOLRayleighPrimary(wavelength,varargin)
p = inputParser;
p.addParameter('p2',false,@(x)(islogical(x)));
p.parse(varargin{:});
p2 = p.Results.p2;

cal = OLGetCalibrationStructure('CalibrationType','BoxBRandomizedLongCableAEyePiece1_12_10_19');

if wavelength == 0 % Make dark spd
    [~,settingsLength] = size(cal.computed.pr650M);
    spd = OLPrimaryToSpd(cal,zeros(settingsLength,1));

else               % Make spd
    % Initial parameters
    fullWidthHalfMax = 20;
    lambda = 0.001;
    
    % Make spd. For p2, divide by 3 
    if p2
        spdRel = OLMakeMonochromaticSpd(cal,wavelength,fullWidthHalfMax)/3;
    else 
        spdRel = OLMakeMonochromaticSpd(cal,wavelength,fullWidthHalfMax);
    end 
    spd = OLMaximizeSpd(cal,spdRel,'lambda',lambda);
end

% Save results
file = sprintf('OLRayleighPrimary_%g.mat',wavelength);
save(fullfile(getpref('ForcedChoiceCM','rayleighDataDir'),...
    'precomputedSpds', file),'spd','p2');
end