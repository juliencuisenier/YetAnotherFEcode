function [C_hcb,Z_hcb] = rotor_forces_reducted(PrimalSub,T_hcb,L_hcb,rho)
%ROTOR_FORCES Summary of this function goes here
%   Omega is the rotationnal speed of the rotor

C_hcb = [];

Z_hcb = [];

for iSub=1:PrimalSub.nSubs
    
    Cs = sparse(PrimalSub.Substructures(iSub).Mesh.nDOFs,PrimalSub.Substructures(iSub).Mesh.nDOFs);

    Zs = sparse(PrimalSub.Substructures(iSub).Mesh.nDOFs,PrimalSub.Substructures(iSub).Mesh.nDOFs);

%     Cs = zeros(PrimalSub.Substructures(iSub).Mesh.nDOFs,PrimalSub.Substructures(iSub).Mesh.nDOFs);
%  
%     Zs = zeros(PrimalSub.Substructures(iSub).Mesh.nDOFs,PrimalSub.Substructures(iSub).Mesh.nDOFs);
    
    for j=1:PrimalSub.Substructures(iSub).Mesh.nNodes
        
        E = sparse(3,PrimalSub.Substructures(iSub).Mesh.nDOFs);
        
        E(:,3*j-2:3*j) = speye(3);
        
        P = [0 1 0; 0 -1 0; 0 0 0];
        
        Cp = rho*P';
        
        Cs = Cs +E'*Cp*E;
        
        J = [1 0 0; 0 1 0; 0  0 0];
        
        Zp = rho*J;
        
        Zs = Zs +E'*Zp*E;
    end
    
    Us_internal = PrimalSub.InternalFreeDOF{iSub};
    Us_interface = PrimalSub.InterfaceDOF{iSub};
    
    Cs_rearranged = [Cs(Us_internal,Us_internal) Cs(Us_internal,Us_interface);...
        Cs(Us_interface,Us_internal) Cs(Us_interface,Us_interface)];
    
    Zs_rearranged = [Zs(Us_internal,Us_internal) Zs(Us_internal,Us_interface);...
        Zs(Us_interface,Us_internal) Zs(Us_interface,Us_interface)];
    
    Cs_hcb = T_hcb{iSub}'*Cs_rearranged*T_hcb{iSub};
    Zs_hcb = T_hcb{iSub}'*Zs_rearranged*T_hcb{iSub};
    
    if isempty(C_hcb)
        C_hcb = L_hcb{iSub}'*Cs_hcb*L_hcb{iSub};
        Z_hcb = L_hcb{iSub}'*Zs_hcb*L_hcb{iSub};
    else
        C_hcb = C_hcb + L_hcb{iSub}'*Cs_hcb*L_hcb{iSub};
        Z_hcb = Z_hcb + L_hcb{iSub}'*Zs_hcb*L_hcb{iSub};
    end



    
end



end

