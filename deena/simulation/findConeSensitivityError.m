function [error,error540] = ...
    findConeSensitivityError(observerParams1,observerParams2,varargin)
% Find root mean square error of the difference in cone spectral
% sensitivities
%
% Syntax:
%   findConeSensitivityError(observerParams1,observerParams2)
%
% Description:
%    Takes two sets of observer parameters and computes the root mean
%    square error of their cone spectral sensitivities, which indicates how
%    similar their perception generally is. Returns two values - a root
%    mean square error across the entire range of cones (typically
%    380-780nm), and a root mean square error at wavelengths greater than
%    540nm (which are typically used in Rayleigh matching). Note that this
%    program ignores observer noise.
%
% Inputs:
%    observerParams1  -Vector of nine individual difference parameters for 
%                      the first observer. See ObserverVecToParams for a 
%                      full description.
%    observerParams1  -Vector of nine individual difference parameters for 
%                      the second observer. See ObserverVecToParams for a 
%                      full description.
%
% Outputs:
%    error       -Root mean square error of the two observers' spectral
%                 sensitivities, computed across the entire range of cones.
%    error540    -Root mean square error of the two observers' spectral
%                 sensitivities, computed at wavelengths greater than
%                 540 nm.
%
% Optional key-value pairs:
%    'S'          -Wavelength sampling for cone calculations, in the 
%                 form [start increment numTerms]. Default is [380 2 201] 
%                 (OneLight convention)  
%    'fieldSize   -Integer field size, in degrees. Default is 2.
%    'age'        -Integer age for simulated observer. Default is 32.
%                        
% History 
%   07/10/20  dce       Wrote it

% Parse input 
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.parse(varargin{:});

% Generate two observers - one with each set of cone parameters 
observer1 = genRayleighObserver('age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'coneVec',observerParams1,'S',p.Results.S);
observer2 = genRayleighObserver('age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'coneVec',observerParams2,'S',p.Results.S);

% Array for storing cone differences. The first column is L cone 
% differences, the second is M differences, and the third is S differences
wls = SToWls(p.Results.S);
coneDiffs = zeros(length(wls),3);

% Calculate squared differences at each wavelength
for i = 1:length(wls)    
    squaredDiffL = (observer1.T_cones(1,i)-...
        observer2.T_cones(1,i))^2;
    squaredDiffM = (observer1.T_cones(2,i)-...
        observer2.T_cones(2,i))^2;
    squaredDiffS = (observer1.T_cones(3,i)-...
        observer2.T_cones(3,i))^2;
    coneDiffs(i,:) = [squaredDiffL squaredDiffM squaredDiffS]; 
end

% Calculate rms error - both overall and at wavelengths over 540 nm 
error = sqrt(mean(coneDiffs,'all'));

ind540 = find(wls==540);
error540 = sqrt(mean(coneDiffs(ind540:end,:),'all')); 
