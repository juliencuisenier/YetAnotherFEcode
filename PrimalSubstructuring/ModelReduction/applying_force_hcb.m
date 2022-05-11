function [F_hcb] = applying_force_hcb(PrimalSub,Fs,T_hcb,L_hcb)
%APPLYING_FORCE_HCB Summary of this function goes here
%   Fs = cell containing forces applied on each substructures

[~,nDOFglobal_hcb] = size(L_hcb{1});

F_hcb = zeros(nDOFglobal_hcb,1);

for iSub=1:PrimalSub.nSubs
   
    
    nInts = length(PrimalSub.InterfaceDOF{iSub});
    nInternalDOF = length(PrimalSub.InternalFreeDOF{iSub});
    
    fs = zeros(nInts + nInternalDOF,1);
    fs(1:nInternalDOF) = Fs{iSub}(PrimalSub.InternalFreeDOF{iSub});
    fs(nInternalDOF+1:end) = Fs{iSub}(PrimalSub.InterfaceDOF{iSub});
    
    fs_hcb = (T_hcb{iSub}'*T_hcb{iSub})\T_hcb{iSub}'*fs;
    F_hcb = F_hcb + L_hcb{iSub}'*fs_hcb;
end

end

