function [Vs_corr] = corr_reduction_indices(PrimalSub,T_hcb,L_hcb,V_hcb)
%CORR_REDUCTIONINDICES Summary of this function goes here
%   Detailed explanation goes here

nSubs = PrimalSub.nSubs;
Vs_corr = cell(1,nSubs);
Vs = cell(1,nSubs);

for iSub=1:nSubs
    Vs_hcb = L_hcb{iSub}*V_hcb;
    Vs{iSub} = T_hcb{iSub}*Vs_hcb;
end

for iSub=1:nSubs
   
    InternalUs = PrimalSub.InternalFreeDOF{iSub};
    InterfaceUs = PrimalSub.InterfaceDOF{iSub};
    Us = PrimalSub.Us{iSub};
    
    nDOFinterface = length(InterfaceUs);
    nDOFinternalFree = length(InternalUs);
    
    
    Vs_corr{iSub} = zeros(length(Us),1);
    
    
    for i=1:nDOFinternalFree
        Vs_corr{iSub}(InternalUs(i)) = Vs{iSub}(i);
    end
    
    for i=1:nDOFinterface
        Vs_corr{iSub}(InterfaceUs(i)) = Vs{iSub}(i+nDOFinternalFree)/2;
    end
    
end

end

