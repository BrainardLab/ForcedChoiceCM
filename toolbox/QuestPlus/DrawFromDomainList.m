function [v] = DrawFromDomainList(domainList)

[vlb, vub] = GetBoundsFromDomainList(domainList);
for ii = 1:length(domainList)
    v(ii) = unifrnd(vlb(ii),vub(ii));   
end