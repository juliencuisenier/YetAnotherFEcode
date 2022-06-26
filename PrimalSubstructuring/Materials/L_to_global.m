function [u] = L_to_global(PrimalSub,us)
%L_TO_GLOBAL Summary of this function goes here
%   From a substructured vector, constructs the global vector, assuming
%   that for a given interface, there are only two substructures

u = zeros(PrimalSub.nDOFglobal,1);

for iSub=1:PrimalSub.nSubs
    us{iSub}(PrimalSub.InterfaceDOF{iSub})=us{iSub}(PrimalSub.InterfaceDOF{iSub})/2;
    u = u + PrimalSub.L{iSub}'*us{iSub};
end

end

