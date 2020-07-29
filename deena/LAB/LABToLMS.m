function [comparisonLMS,referenceLMS,comparisonXYZ,referenceXYZ]...
    = LABToLMS(comparisonLab,referenceLab,referenceXYZ,varargin)
% Computes cone responses from CIELAB coordinates
% Syntax:
%   LABToLMS(comparisonLab,referenceLab,referenceXYZ)
%
% Description:
%    Takes in a pair of CIELAB coordinates and a reference XYZ triplet. 
%    Converts the CIELAB coordinates first to XYZ and then to LMS cone 
%    responses. 
%
% Inputs:
%    comparisonLab   -3x1 CIELAB coordinate of comparison light
%    referenceLab    -3x1 CIELAB coordinate of reference light 
%    referenceXYZ    -3x1 XYZ coordinate of reference light 
%
% Outputs:
%    comparisonLMS   -3x1 LMS response for the first CIELAB coordinate 
%    referenceLMS    -3x1 LMS response for the second CIELAB coordinate 
%    comparisonXYZ   -3x1 XYZ coordinates for the first CIELAB coordinate 
%    referenceXYZ    -3x1 XYZ coordinates for the second CIELAB coordinate 
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
%    dce    7/29/20   -Wrote it 

% Parse input 
p = inputParser;
p.addParameter('T_cones',[],@(x)(isnumeric(x)));
p.addParameter('S_cones',[],@(x)(isnumeric(x)));
p.addParameter('T_xyz',[],@(x)(isnumeric(x)));
p.addParameter('S_xyz',[],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Lab to XYZ
referenceXYZ = LabToXYZ(referenceLab,referenceXYZ);
comparisonXYZ = LabToXYZ(comparisonLab,referenceXYZ);

% XYZ to LMS 
% Load the appropriate cone parameters
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

% Conversion matrices
M_LMSToXYZ = (T_cones'\T_xyz')';
M_XYZToLMS = inv(M_LMSToXYZ);

% LMS values
comparisonLMS = M_XYZToLMS*comparisonXYZ;    
referenceLMS = M_XYZToLMS*referenceXYZ;  
end
