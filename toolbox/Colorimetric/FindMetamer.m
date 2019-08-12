function [metamerPrimary,metamerSpd,metamerLMS] = FindMetamer(apparatusParams,T,spectrumLMS)
% Find apparatus metamer for a passed spectrum with respect to passed fundamentals.
%
% Syntax:
%    [metamerPrimary,metamerSpd,metamerLMS] = FindMetamer(apparatusParams,T,spectrumLMS)   
%
% Description:
%    Find a metamer for passed spectrum within apparatus capability.
%
% Inputs:
%    apparatusParams                       - Structure describing apparatus.    
%    T                                     - Cone fundametals
%    spectrumLMS                           - LMS of spectrum for which to find metamer.
%
% Outputs:
%    metamerPrimary                        - Weights on appartus primaries for metamer.
%    metamerSpd                            - Spectrum of metamer
%    metamerLMS                            - Metamer's LMS cone coordinates.
%
% Optional key/value pairs:
%
% See also: DefaultApparatusParams, LMSToPrimary
%

% History:
%   08/11/19  dhb  Wrote it.

switch (apparatusParams.type)
    case 'monochromatic'
        
        % Find apparatus spectrum that hits it
        [metamerPrimary,metamerSpd] = LMSToPrimary(apparatusParams,T,spectrumLMS);

        % Check
        metamerLMS = T*metamerSpd;
        if (max(abs(metamerLMS-spectrumLMS)./spectrumLMS) > 1e-6)
            error('Failed to compute good metamer');
        end
        
    otherwise
        error('Unknown apparatus type specified in parameter structure');
end
