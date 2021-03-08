function [refIntensity, primaryRatio] = OLSpdToPittPoint(refSpd,primarySpd,darkSpd,lightFileName)
% Short program for computing primary ratio and reference intensity from
% match spds
%
% Syntax:
%   OLSpdToPittPoint(refSpd,primarySpd,darkSpd,lightFile)
%
% Description:
%    Takes in a pair of spds that were identified as a match. Using 
%    regression, computes the primary ratio and test intensity relative to 
%    the nominal spectra at full power.
%
% Inputs:
%    refSpd           - Reference spd, entered as a column vector.
%    primarySpd       - Primary mixture spd, entered as a column vector.
%    darkSpd          - Dark spd, entered as a column vector.
%    lightFileName    - Full filepath to the data file containing baseline
%                       nominal spds.
% Outputs:
%    refIntensity     - Estimated reference intensity, between 0 and 1.
%    primaryRatio     - Estimated primary ratio, between 0 and 1.
%
% Optional key-value pairs:
%    None

% History:
%   2/25/21  dce       Wrote it

% Subtract dark spds from the measured lights
refSpdND = refSpd - darkSpd;
primarySpdND = primarySpd - darkSpd;

% Nominal spds for the lights at their peak power. These have had dark spds
% subtracted. 
lightFile = load(lightFileName);
nominalP1Spd = lightFile.primary1IncrSpdPredicted*lightFile.p1ScaleFactor;
nominalP2Spd = lightFile.primary2IncrSpdPredicted*lightFile.p2ScaleFactor;
nominalRefSpd = lightFile.testIncrSpdPredicted*lightFile.testScaleFactor;

% Estimate reference intensity using regression 
refIntensity = nominalRefSpd\refSpdND;

% Estimate primary ratio using regression 
%
% We want primarySpdND = nominalP1Spd*primaryRatio + nominalP2Spd*(1-primaryRatio)
% This leads to the following  equation: 
% primarySpdND = nominalP1Spd*primaryRatio - nominalP2Spd*primaryRatio + nominalP2Spd
% [primarySpdND-nominalP2Spd] = primaryRatio[nominalP1Spd-nominalP2Spd]
% primaryRatio = [nominalP1Spd-nominalP2Spd]\[primarySpdND-nominalP2Spd]

primaryRatio = (nominalP1Spd-nominalP2Spd)\(primarySpdND-nominalP2Spd);

% Check 
makePlots = false;
if makePlots
    wls = SToWls([380 2 201]);
    figure();
    subplot(2,1,1);
    hold on;
    plot(wls,primarySpdND,'LineWidth',3);
    plot(wls,nominalP1Spd*primaryRatio+nominalP2Spd*(1-primaryRatio),'LineWidth',1.5);
    legend('Measured Spd', 'Calculated Spd');
    title('Primary Mixture');
    xlabel('Wavelength (nm)');
    ylabel('Power');
    
    subplot(2,1,2);
    hold on;
    plot(wls,refSpdND,'LineWidth',3);
    plot(wls,nominalRefSpd*refIntensity,'LineWidth',1.5);
    legend('Measured Spd', 'Fitted Spd');
    title('Reference Light')
    xlabel('Wavelength (nm)');
    ylabel('Power');
end 

end 