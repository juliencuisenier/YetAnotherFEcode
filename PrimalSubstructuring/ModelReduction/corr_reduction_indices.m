function [V_corr] = corr_reduction_indices(PrimalSub,V)
%CORR_REDUCTIONINDICES Summary of this function goes here
%   Detailed explanation goes here

V_corr = zeros(length(V),1);
InternalDOFDone = 0;

for iSub=1:PrimalSub.nSubs
    
    V_corr(PrimalSub.InternalFreeDOF{iSub}) = V(InternalDOFDone+1:InternalDOFDone+length(PrimalSub.InternalFreeDOF{iSub}));
    InternalDOFDone = InternalDOFDone + length(PrimalSub.InternalFreeDOF{iSub});
    
end

InterfaceDOFDone = 0;

for iSub=1:PrimalSub.nSubs
    
    V_corr(PrimalSub.InterfaceDOF{iSub}) = V(InternalDOFDone+InterfaceDOFDone+1:InternalDOFDone+InterfaceDOFDone+length(PrimalSub.InterfaceDOF{iSub}));
    InterfaceDOFDone = InterfaceDOFDone + length(PrimalSub.InterfaceDOF{iSub});
    
end

end

