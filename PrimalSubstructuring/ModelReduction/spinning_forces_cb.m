function [G_cb,Ksp_cb] = spinning_forces_cb(PrimalSub,T_cb,L_cb,Omega)
%ROTOR_FORCES_CB Summary of this function goes here
%   Compute the coriolis and spin softening matrices for a Craig-Bampton
%   reduced model

G_cb = [];
Ksp_cb = [];

for iSub=1:PrimalSub.nSubs
    Us_internal = PrimalSub.InternalFreeDOF{iSub};
    Us_interface = PrimalSub.InterfaceDOF{iSub};
    
    iG = PrimalSub.Substructures(iSub).coriolis_matrix(Omega);
    iKsp = PrimalSub.Substructures(iSub).spin_softening_matrix(Omega);
    
    iG = [iG(Us_internal,Us_internal) iG(Us_internal,Us_interface);...
        iG(Us_interface,Us_internal) iG(Us_interface,Us_interface)];
    
    iKsp = [iKsp(Us_internal,Us_internal) iKsp(Us_internal,Us_interface);...
        iKsp(Us_interface,Us_internal) iKsp(Us_interface,Us_interface)];
    
    iG_cb = T_cb{iSub}'*iG*T_cb{iSub};
    iKsp_cb = T_cb{iSub}'*iKsp*T_cb{iSub};
    
    if isempty(G_cb)
        G_cb = L_cb{iSub}'*iG_cb*L_cb{iSub};
        Ksp_cb = L_cb{iSub}'*iKsp_cb*L_cb{iSub};
    else
        G_cb = G_cb + L_cb{iSub}'*iG_cb*L_cb{iSub};
        Ksp_cb = Ksp_cb + L_cb{iSub}'*iKsp_cb*L_cb{iSub};
    end
end
end

