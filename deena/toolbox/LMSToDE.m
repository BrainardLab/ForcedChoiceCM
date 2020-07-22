function [comparisonDE,referenceLab,comparisonLab] = ...
    LMSToDE(referenceLMS,comparisonLMS,varargin)
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
%    referenceLMS    -3x1 cone responses for first light 
%    comparisonLMS   -3x1 cone responses for second light 
%
% Outputs:
%    comparisonDE    -Scalar value of delta E for the two lights.
%    referenceLab    -3x1 CIELAB coordinates for the first spectrum.
%    comparisonLab   -3x1 CIELAB coordinates for the second spectrum.
%
% Optional key-value pairs:
%    'T_cones'       -Passed cone fundamentals, in a 3xn matrix. Default is 
%                     [], in which case the fundamentals are loaded from 
%                     'T_cones_ss2'.
%    'S_cones'       -Wavelength sampling for the cones, in the form 
%                     [start delta nterms]. Default is []. Must be
%                     specified if T_cones is nonempty.
%    'T_xyz'         -Passed xyz fundamentals, in a 3xm matrix. Default 
%                     is [], in which case the fundamentals are loaded from 
%                     'T_xyz1931'.
%    'S_xyz'         -Wavelength sampling for the cones, in the form 
%                     [start delta nterms]. Default is []. Must be
%                     specified if T_xyz is nonempty.

% History
%    dce    7/20/20   -Adapted from example code from dhb
%    dce    7/22/20   -Added options to specify cone input

% Parse input 
p = inputParser;
p.addParameter('T_cones',[],@(x)(isnumeric(x)));
p.addParameter('S_cones',[],@(x)(isnumeric(x)));
p.addParameter('T_xyz',[],@(x)(isnumeric(x)));
p.addParameter('S_xyz',[],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Load cone data 
if isempty(p.Results.T_cones)
    cones = load('T_cones_ss2');
    T_cones = cones.T_cones_ss2;
    S_cones = cones.S_cones_ss2;
else 
    if isempty(p.Results.S_cones)
        error('S must be provided along with cone fundamentals'); 
    else 
        T_cones = p.Results.T_cones;
        S_cones = p.Results.S_cones;
    end 
end 

% Load XYZ data
if isempty(p.Results.T_xyz)
    xyzData = load('T_xyz1931');
    T_xyz = SplineCmf(xyzData.S_xyz1931,xyzData.T_xyz1931,S_cones);
else 
    if isempty(p.Results.S_xyz)
        error('S must be provided along with xyz functions'); 
    else 
        T_xyz = SplineCmf(p.Results.S_xyz,p.Results.T_xyz,S_cones);
    end 
end 
M_LMSToXYZ = (T_cones'\T_xyz')';

% Check if the matrix properly recovers XYZ functions from cone functions
CHECK_MATRIX = false;
if (CHECK_MATRIX)
    T_xyzCheck = M_LMSToXYZ*T_cones;
    figure; clf; hold on
    plot(T_xyz','r','LineWidth',4);
    plot(T_xyzCheck','b','LineWidth',2);
end

% Compute XYZ coordinates
referenceXYZ = M_LMSToXYZ*referenceLMS;
comparisonXYZ = M_LMSToXYZ*comparisonLMS;

% Compute CIELAB coordinates and difference
referenceLab = XYZToLab(referenceXYZ,referenceXYZ);
comparisonLab = XYZToLab(comparisonXYZ,referenceXYZ);
comparisonDE = ComputeDE(comparisonLab,referenceLab);
end