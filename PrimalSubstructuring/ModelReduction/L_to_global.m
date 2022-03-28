function [u] = L_to_global(PrimalSub,us)
%L_TO_GLOBAL uses the localization matrix to have a global vector from a
%vector defined on each substructure

%   PrimalSub : PrimalSubstructure class
%   us : a cell with a length of nSubs

u = [];

for iSub=1:PrimalSub.nSubs
    if isempty(u)
        u = PrimalSub.L{iSub}'*us{iSub};
        
    else
        u = u + PrimalSub.L{iSub}'*us{iSub};
    end
end

end

