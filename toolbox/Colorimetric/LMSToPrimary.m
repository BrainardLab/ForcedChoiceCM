function [spectrumPrimary,spectrum] = LMSToPrimary(apparatusParams,T,spectrumLMS)
% Get primary weights and spectrum that have passed LMS coordinates
%
% Syntax:
%    [spectrumPrimary,spectrum] = LMSToPrimary(apparatusParams,T,spectrumLMS)
%
% Description:
%    Get apparatus primary weights to produce spectrum with desired LMS
%    coordinates.  Also returns the spectrum.
%
% Inputs:
%    apparatusParams                       - Structure describing apparatus.    
%    T                                     - Cone fundametals
%    spectrumLMS                           - Target LMS coordinates.
%
% Outputs:
%    spectrumPrimary                       - Spectrum primary weights.
%    spectrum                              - Spectrum of metamer%
%
% Optional key/value pairs:
%    None.
%
% See also: DefaultApparatusParams, FindMetamer, PrimaryToLMS
%

% History:
%   08/11/19  dhb  Wrote it.

switch (apparatusParams.type)
    case 'monochromatic'
        
        % Get conversion matrices
        M_PrimaryToLMS = T*apparatusParams.primaryBasis;
        M_LMSToPrimary = inv(M_PrimaryToLMS);
        
        % Convert
        spectrumPrimary = M_LMSToPrimary*spectrumLMS;
        spectrum = apparatusParams.primaryBasis*spectrumPrimary;
        checkLMS = T*spectrum;
        if (max(abs(checkLMS-spectrumLMS)./spectrumLMS) > 1e-6)
            error('Failed to reproduce desired LMS');
        end
        
    otherwise
        error('Unknown apparatus type specified in parameter structure');
end


