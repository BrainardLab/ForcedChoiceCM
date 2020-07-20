function [comparisonDE,referenceLab,comparisonLab] = LMSToDE(referenceLMS,comparisonLMS,varargin)
% Computes CIELAB delta E based on cone responses to two spectra

% Syntax:
%   LMSToDE(referenceLMS,comparisonLMS)
%
% Description:
%    Takes in cone responses to two spectra, and converts these first to
%    XYZ coordinates and then to CIELAB coordinates. Returns delta E, a
%    measure of how far apart the two CIELAB coordinates are. 
%
% Inputs:
%    referenceLMS    -LMS response to the first spectrum 
%    comparisonLMS   -LMS response to the second spectrum
%
% Outputs:
%    comparisonDE    -Scalar value of delta E for the two lights.
%    referenceLab    -CIELAB coordinates for the first light.
%    comparisonLab   -CIELAB coordinates for the second light.
%
% Optional key-value pairs:
%    'coneDataFile'    -Character vector of filename for cone data. Default
%                       is 'T_cones_ss2'. 
%    'xyzDataFile'     -Character vector of filename for xyz function data.
%                       Default is 'T_xyz1931'.

% History
%    dce    7/29/20   -Adapted from example code from dhb

p = inputParser;
p.addParameter('coneDataFile','T_cones_ss2',@(x)(ischar(x)));
p.addParameter('xyzDataFile','T_xyz1931',@(x)(ischar(x)));
p.parse(varargin{:});

cones = load(coneDataFile);
T_cones = cones.T_cones_ss2;
S = S_cones_ss2;

xyzData = load(xyzDataFile);
T_xyz = SplineCmf(xyzData.S_xyz1931,xyzData.T_xyz1931,S);
M_LMSToXYZ = (T_cones'\T_xyz')';

CHECK_MATRIX = false;
if (CHECK_MATRIX)
    T_xyzCheck = M_LMSToXYZ*T_cones;
    figure; clf; hold on
    plot(T_xyz','r','LineWidth',4);
    plot(T_xyzCheck','b','LineWidth',2);
end

referenceXYZ = M_LMSToXYZ*referenceLMS;
comparisonXYZ = M_LMSToXYZ*comparisonLMS;

referenceLab = XYZToLab(referenceXYZ,referenceXYZ);
comparisonLab = XYZToLab(comparisonXYZ,referenceXYZ);
comparisonDE = ComputeDE(comparisonLab,referenceLab);
end