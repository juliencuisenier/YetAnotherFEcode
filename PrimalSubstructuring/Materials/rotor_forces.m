function [Cr,Z] = rotor_forces(PrimalSub,rho)
%ROTOR_FORCES Summary of this function goes here
%   Omega is the rotationnal speed of the rotor

Cr = sparse(PrimalSub.nDOFglobal,PrimalSub.nDOFglobal);
Z = sparse(PrimalSub.nDOFglobal,PrimalSub.nDOFglobal);

Zs = cell(1,PrimalSub.nSubs);
Cs = cell(1,PrimalSub.nSubs);

for iSub=1:PrimalSub.nSubs
    
    C_s = sparse(PrimalSub.Substructures(iSub).Mesh.nDOFs,PrimalSub.Substructures(iSub).Mesh.nDOFs);

    Z_s = sparse(PrimalSub.Substructures(iSub).Mesh.nDOFs,PrimalSub.Substructures(iSub).Mesh.nDOFs);
    
    for j=1:PrimalSub.Substructures(iSub).Mesh.nNodes
        
        E = sparse(3,PrimalSub.Substructures(iSub).Mesh.nDOFs);
        
        E(:,3*j-2:3*j) = speye(3);
        
        P = [0 1 0; 0 -1 0; 0 0 0];
        
        Cp = rho*P';
        
        C_s = C_s +E'*Cp*E;
        
        J = [1 0 0; 0 1 0; 0  0 0];
        
        Zp = rho*J;
        
        Z_s = Z_s +E'*Zp*E;
    end
    
    Cr = Cr + PrimalSub.L{iSub}'*C_s*PrimalSub.L{iSub};  
    Z = Z + PrimalSub.L{iSub}'*Z_s*PrimalSub.L{iSub};
    
end


end

