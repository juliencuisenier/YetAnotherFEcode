function [M_cb,C_cb,K_cb,T_cb,L_cb] = CraigBamptonReduction(PrimalSub,Omega,freq)
%GREGBAMPTONREDUCTION Summary of this function goes here
%   Detailed explanation goes here

nSubs = PrimalSub.nSubs;

T_cb = cell(1,nSubs);

Us_cb = cell(1,nSubs);
L_cb = cell(1,nSubs);

Mcs_rearranged = cell(1,nSubs);
Kcs_rearranged = cell(1,nSubs);
Gcs_rearranged = cell(1,nSubs);

ms = zeros(1,nSubs);

nInternalDOFglobal_cb = 0;

nDOFPerNode = PrimalSub.Substructures(1).Mesh.nDOFPerNode;

%First loop to determine the T_cb matrix
for iSub=1:nSubs
    
    Us_internal = PrimalSub.InternalFreeDOF{iSub};
    Us_interface = PrimalSub.InterfaceDOF{iSub};

    nInts = length(Us_interface);
    
%% Coriolis matrix

% These matrices will be used later, if Omega isn't empty


%% Component modes
    u0 = zeros(length(PrimalSub.Us{iSub}),1);
    [Ks,~] = PrimalSub.Substructures(iSub).tangent_stiffness_and_force(u0);
    
    % Spinning force matrices, if the structure is rotated
    if ~isempty(Omega)
        Gs = PrimalSub.Substructures(iSub).coriolis_matrix(Omega);
        Gcs_rearranged{iSub} = [Gs(Us_internal,Us_internal) Gs(Us_internal,Us_interface);...
            Gs(Us_interface,Us_internal) Gs(Us_interface,Us_interface)];
        
        Ks = Ks + PrimalSub.Substructures(iSub).spin_softening_matrix(Omega);        
    end
    
    Ks_ib = Ks(Us_internal,Us_interface);
    Ks_ii = Ks(Us_internal,Us_internal);
    CM = [-Ks_ii\Ks_ib;eye(nInts)];
    
    
%% Internal VM
    Ms = PrimalSub.Substructures(iSub).mass_matrix();
    Ms_ii = Ms(Us_internal,Us_internal);
    
    [V0_i,om] = eigs(Ks_ii, Ms_ii, 100, 'SM'); 
    f = sqrt(diag(om))/2/pi;
    f = nonzeros(f < 2*freq);
    [m,~] = size(f);
    [~,ind] = sort(f);
    V0_i = V0_i(:,ind);
    
    IVM = [V0_i;zeros(nInts,m)];
    
    T_cb{iSub} = [IVM CM]; 
    Mcs_rearranged{iSub} = [Ms_ii Ms(Us_internal,Us_interface); Ms(Us_interface,Us_internal) Ms(Us_interface,Us_interface)];
    Kcs_rearranged{iSub} = [Ks_ii Ks_ib; Ks(Us_interface,Us_internal) Ks(Us_interface,Us_interface)];
    ms(iSub) = m;
    
    Us_cb{iSub} = zeros(m+nInts,1);
    
    nInternalDOFglobal_cb = nInternalDOFglobal_cb + m;
end

nDOFglobal_cb = nInternalDOFglobal_cb + PrimalSub.nInt*nDOFPerNode;

U = 1:nDOFglobal_cb';


InterfaceDOFdone = nInternalDOFglobal_cb;
InternalDOFdone = 0;
IntsDOFdone = zeros(1,nSubs);


Interfaces = PrimalSub.Interfaces;

%Second Loop to assignate the DOFs
for iSub=1:nSubs

    m = ms(iSub);
    Us_cb{iSub}(1:m) = U(InternalDOFdone+1:InternalDOFdone+m);
    
    InternalDOFdone = InternalDOFdone + m;
    
    Interfaces_loc = Interfaces(Interfaces(:,iSub)~=0,:);
    
    for jSub=iSub+1:nSubs
        ind = find(Interfaces_loc(:,jSub));
        if ~isempty(ind)
            nInts_loc = length(ind)*nDOFPerNode;
            Us_cb{iSub}(m+IntsDOFdone(iSub)+1:m+IntsDOFdone(iSub)+nInts_loc) = U(InterfaceDOFdone+1:InterfaceDOFdone+nInts_loc);
            Us_cb{jSub}(ms(jSub)+IntsDOFdone(jSub)+1:ms(jSub)+IntsDOFdone(jSub)+nInts_loc) = U(InterfaceDOFdone+1:InterfaceDOFdone+nInts_loc);
            InterfaceDOFdone = InterfaceDOFdone + nInts_loc;
            IntsDOFdone(iSub) = IntsDOFdone(iSub) + nInts_loc;
            IntsDOFdone(jSub) = IntsDOFdone(jSub) + nInts_loc;
        end
    end 
    
    
    
end

%Loop to build L_cb

for iSub=1:nSubs  
    us = Us_cb{iSub}; %DOF vector of substructure i
    n_s = length(us);
    
    Ls = sparse(1:n_s, us, true(n_s,1), n_s, nDOFglobal_cb );
    L_cb{iSub} = Ls;
end

M_cb = [];
K_cb = [];

for iSub=1:nSubs
      
    Kcs = Kcs_rearranged{iSub};
    Mcs = Mcs_rearranged{iSub};
    Gcs = Gcs_rearranged{iSub};
    
    
    Ms_cb = T_cb{iSub}'*Mcs*T_cb{iSub};
    Gs_cb = T_cb{iSub}'*Gcs*T_cb{iSub};
    Ks_cb = T_cb{iSub}'*Kcs*T_cb{iSub};
    
    if isempty(M_cb)
        M_cb = L_cb{iSub}'*Ms_cb*L_cb{iSub};
        C_cb = L_cb{iSub}'*Gs_cb*L_cb{iSub};
        K_cb = L_cb{iSub}'*Ks_cb*L_cb{iSub};
    else
        M_cb = M_cb + L_cb{iSub}'*Ms_cb*L_cb{iSub};
        C_cb = C_cb + L_cb{iSub}'*Gs_cb*L_cb{iSub};
        K_cb = K_cb + L_cb{iSub}'*Ks_cb*L_cb{iSub};
    end
        
end


