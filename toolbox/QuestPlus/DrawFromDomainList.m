function [v] = DrawFromDomainList(domainList)
% Draw parameters from cell array of domain for each parameter
%
% Syntax:
%   [v] = DrawFromDomainList(domainList)
%
% Description:
%    Determine lower and upper parameter bound for each parameter from the
%    domain list, and then draw uniformly from the range.
%
% Inputs:
%    domainList                           - Cell array where each entry is
%                                           the domain for the
%                                           corresponding parameter, in the
%                                           form used by qpInitialize.
%
% Outputs:
%    v                                    - The random draw in row vector
%                                           form.
%
% Optional key/value pairs:
%    None.
%
% See also: GetBoundsFromDomainList, qpInitialize
%

% History:
%   08/25/19  dhb  Wrote it

% Examples:
%{
    psiParamsDomainList = {0, 0, -20:5:20, -20:5:20, 0, -4:1:4, -4:1:4, 0, 0.02};
    v = DrawFromDomainList(psiParamsDomainList)
%}

[vlb, vub] = GetBoundsFromDomainList(domainList);
for ii = 1:length(domainList)
    v(ii) = unifrnd(vlb(ii),vub(ii));   
end