function [error,LError,MError,SError] = ...
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
%    similar their perception generally is. Returns a root mean square 
%    error across the entire range of cones, as well as individual error
%    terms for L, M, and S. 
%
% Inputs:
%    observerParams1  -Vector of eight individual difference parameters for 
%                      the first observer. See ObserverVecToParams for a 
%                      full description.
%    observerParams1  -Vector of eight individual difference parameters for 
%                      the second observer. See ObserverVecToParams for a 
%                      full description.
%
% Outputs:
%    error       -Root mean square error of the two observers' spectral
%                 sensitivities, computed across the entire range of cones.
%    LError       -Root mean square error of the two observers' L cone 
%                 spectral sensitivities.
%    MError      -Root mean square error of the two observers' M cone 
%                 spectral sensitivities.
%    SError      -Root mean square error of the two observers' S cone 
%                 spectral sensitivities. 
%
% Optional key-value pairs:
%    'S'          -Wavelength sampling for cone calculations, in the 
%                 form [start increment numTerms]. Default is [380 2 201] 
%                 (OneLight convention)  
%    'fieldSize'  -Integer field size, in degrees. Default is 2.
%    'age'        -Integer age for simulated observer. Default is 32.
%    'opponentParams' -4-element vector with opponent contrast parameters. 
%                      (1) is the luminance weight, (2) is the rg weight,
%                      (3) is the by weight, and (4) is the noise standard
%                      deviation. Default is [40.3908 205.7353 62.9590 1.0000].
%                        
% History 
%   07/10/20  dce       Wrote it
%   08/05/20  dce       Added opponent cone parameters.
%   12/22/20  dce       Added individual cone errors as additional outputs.

% Parse input 
p = inputParser;
p.addParameter('age',32,@(x)(isnumeric(x)));
p.addParameter('fieldSize',2,@(x)(isnumeric(x)));
p.addParameter('S',[380 2 201],@(x)(isnumeric(x)));
p.addParameter('opponentParams',[40.3908 205.7353 62.9590 1.0000],@(x)(isvector(x)));
p.parse(varargin{:});

% Generate two observers - one with each set of cone parameters 
observer1 = genRayleighObserver('age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'coneVec',observerParams1,'S',p.Results.S,...
    'opponentParams',p.Results.opponentParams);
observer2 = genRayleighObserver('age',p.Results.age,'fieldSize',...
    p.Results.fieldSize,'coneVec',observerParams2,'S',p.Results.S,...
    'opponentParams',p.Results.opponentParams);

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

% Calculate rms error 
errorChannels = sqrt(mean(coneDiffs,1));
LError = errorChannels(1);
MError = errorChannels(2);
SError = errorChannels(3);
error = sqrt(mean(coneDiffs,'all'));