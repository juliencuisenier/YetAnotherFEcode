function [F_cb] = applying_force_cb(PrimalSub,Fs,T_cb,L_cb)
%APPLYING_FORCE_HCB Summary of this function goes here
%   Fs = cell containing forces applied on each substructures

nDOFglobal_hcb = size(L_cb{1},2);

F_cb = zeros(nDOFglobal_hcb,1);

for iSub=1:PrimalSub.nSubs
     
    nInts = length(PrimalSub.InterfaceDOF{iSub});
    nInternalDOF = length(PrimalSub.InternalFreeDOF{iSub});
    
    fs = zeros(nInts + nInternalDOF,1);
    fs(1:nInternalDOF) = Fs{iSub}(PrimalSub.InternalFreeDOF{iSub});
    fs(nInternalDOF+1:end) = Fs{iSub}(PrimalSub.InterfaceDOF{iSub});
    
    fs_cb = T_cb{iSub}'*fs;
    F_cb = F_cb + L_cb{iSub}'*fs_cb;
end

end

