function [us] = L_to_local(PrimalSub,u)
%L_TO_LOCAL uses the localization matrix to have a vector on each substructure
%from a global vector

%   PrimalSub : PrimalSubstructure class
%   u : a vector

us = {};

for iSub=1:PrimalSub.nSubs
    us{iSub} =PrimalSub.L{iSub}*u;
end

end

