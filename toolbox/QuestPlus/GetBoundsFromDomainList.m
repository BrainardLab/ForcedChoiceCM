function [vlb, vub] = GetBoundsFromDomainList(domainList)

for ii = 1:length(domainList)
    vlb(ii) = min(domainList{ii});
    vub(ii) = max(domainList{ii});
end