function theFig = OLPlotSpdCheck(wls, spds)
% Plots an array of target spds 
% Syntax:
%   OLPlotSpdCheck
%
% Description:
%    Plots target spds across a given range of wavelengths, producing a
%    single figure
%
% Inputs:
%    wls    - nx1 array of wavelengths (nm)
%    spds   - nxm array of spds (each one over the wavelengths in wls)
%
% Outputs:
%    theFig - Figure handle to plot 
%
% Optional key-value pairs:
%    none 

% History:
%   de  11/26/19  - modified script
%   de  3/29/20   - added documentation and figure handle output
%   de  7/16/20   - improved formatting (line thickness)

% wls = 380:2:780; 
[~, col] = size(spds);
theFig = figure(); 
lineWidthBase = 2; 
lineWidthInc = lineWidthBase/col; 
hold on; 
for i = 1:col
    plot(wls,spds(:,i),'LineWidth',lineWidthBase-lineWidthInc*(i-1));
end
xlabel('Wavelength(nm)');
ylabel('Power');
hold off; 
end 