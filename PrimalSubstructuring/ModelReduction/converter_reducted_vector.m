function [V_corr] = converter_reducted_vector(PrimalSub,T_cb,L_cb,V_cb)
%CORR_REDUCTIONINDICES Summary of this function goes here
%   Detailed explanation goes here

nSubs = PrimalSub.nSubs;
Vs = cell(1,nSubs);
Vs_corr = cell(1,nSubs);

for iSub=1:nSubs
    Vs_hcb = L_cb{iSub}*V_cb;
    Vs{iSub} = T_cb{iSub}*Vs_hcb;
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
        Vs_corr{iSub}(InterfaceUs(i)) = Vs{iSub}(i+nDOFinternalFree);
    end  
end

V_corr = L_to_global(PrimalSub,Vs_corr);

end

