function [M_hcb,K_hcb,T_hcb,L_hcb] = CraigBamptonReduction(PrimalSub,freq)
%GREGBAMPTONREDUCTION Summary of this function goes here
%   Detailed explanation goes here

nSubs = PrimalSub.nSubs;

T_hcb = cell(1,nSubs);

Us_hcb = cell(1,nSubs);
L_hcb = cell(1,nSubs);

Mcs_rearranged = cell(1,nSubs);
Kcs_rearranged = cell(1,nSubs);

ms = zeros(1,nSubs);

nInternalDOFglobal_hcb = 0;

nDOFPerNode = PrimalSub.Substructures(1).Mesh.nDOFPerNode;

%First loop to determine the T_hcb matrix
for iSub=1:nSubs
    
    Us_internal = PrimalSub.InternalFreeDOF{iSub};
    Us_interface = PrimalSub.InterfaceDOF{iSub};

    nInts = length(Us_interface);
%% Component modes
    u0 = zeros(length(PrimalSub.Us{iSub}),1);
    [Ks,~] = PrimalSub.Substructures(iSub).tangent_stiffness_and_force(u0);
    Ks_ib = Ks(Us_internal,Us_interface);
    Ks_ii = Ks(Us_internal,Us_internal);
    CM = [-Ks_ii\Ks_ib;eye(nInts)];
    
    
%% Internal VM
    Ms = PrimalSub.Substructures(iSub).mass_matrix();
    Ms_ii = Ms(Us_internal,Us_internal);
    
    [V0_i,om] = eigs(Ks_ii, Ms_ii, 100, 'SM'); 
    f = sqrt(diag(om))/2/pi;
    f = nonzeros(f < freq);
    [m,~] = size(f);
    [~,ind] = sort(f);
    V0_i = V0_i(:,ind);
    
    IVM = [V0_i;zeros(nInts,m)];
    
    T_hcb{iSub} = [IVM CM]; 
    Mcs_rearranged{iSub} = [Ms_ii Ms(Us_internal,Us_interface); Ms(Us_interface,Us_internal) Ms(Us_interface,Us_interface)];
    Kcs_rearranged{iSub} = [Ks_ii Ks_ib; Ks(Us_interface,Us_internal) Ks(Us_interface,Us_interface)];
    ms(iSub) = m;
    
    Us_hcb{iSub} = zeros(m+nInts,1);
    
    nInternalDOFglobal_hcb = nInternalDOFglobal_hcb + m;
end

nDOFglobal_hcb = nInternalDOFglobal_hcb + PrimalSub.nInt*nDOFPerNode;

U = 1:nDOFglobal_hcb';


InterfaceDOFdone = nInternalDOFglobal_hcb;
InternalDOFdone = 0;
IntsDOFdone = zeros(1,nSubs);


Interfaces = PrimalSub.Interfaces;

%Second Loop to assignate the DOFs
for iSub=1:nSubs

    m = ms(iSub);
    Us_hcb{iSub}(1:m) = U(InternalDOFdone+1:InternalDOFdone+m);
    
    InternalDOFdone = InternalDOFdone + m;
    
    Interfaces_loc = Interfaces(Interfaces(:,iSub)~=0,:);
    
    for jSub=iSub+1:nSubs
        ind = find(Interfaces_loc(:,jSub));
        if ~isempty(ind)
            nInts_loc = length(ind)*nDOFPerNode;
            Us_hcb{iSub}(m+IntsDOFdone(iSub)+1:m+IntsDOFdone(iSub)+nInts_loc) = U(InterfaceDOFdone+1:InterfaceDOFdone+nInts_loc);
            Us_hcb{jSub}(ms(jSub)+IntsDOFdone(jSub)+1:ms(jSub)+IntsDOFdone(jSub)+nInts_loc) = U(InterfaceDOFdone+1:InterfaceDOFdone+nInts_loc);
            InterfaceDOFdone = InterfaceDOFdone + nInts_loc;
            IntsDOFdone(iSub) = IntsDOFdone(iSub) + nInts_loc;
            IntsDOFdone(jSub) = IntsDOFdone(jSub) + nInts_loc;
        end
    end 
    
    
    
end

%Loop to build L_hcb

for iSub=1:nSubs  
    us = Us_hcb{iSub}; %DOF vector of substructure i
    n_s = length(us);
    
    Ls = sparse(1:n_s, us, true(n_s,1), n_s, nDOFglobal_hcb );
    L_hcb{iSub} = Ls;
end

M_hcb = [];
K_hcb = [];

for iSub=1:nSubs
      
    Kcs = Kcs_rearranged{iSub};
    Mcs = Mcs_rearranged{iSub};
    
    
    Ms_hcb = T_hcb{iSub}'*Mcs*T_hcb{iSub};
    
    Ks_hcb = T_hcb{iSub}'*Kcs*T_hcb{iSub};
    
    if isempty(M_hcb)
        M_hcb = L_hcb{iSub}'*Ms_hcb*L_hcb{iSub};
        K_hcb = L_hcb{iSub}'*Ks_hcb*L_hcb{iSub};
    else
        M_hcb = M_hcb + L_hcb{iSub}'*Ms_hcb*L_hcb{iSub};
        K_hcb = K_hcb + L_hcb{iSub}'*Ks_hcb*L_hcb{iSub};
    end
        
end

end


