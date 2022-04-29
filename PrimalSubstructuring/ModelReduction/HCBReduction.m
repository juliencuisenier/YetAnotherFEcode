function [M_hcb,K_hcb,L_hcb] = HCBReduction(PrimalSub,x,nM)
%HCBREDUCTION Summary of this function goes here

nSubs = PrimalSub.nSubs;

Us_hcb = cell(1,nSubs);      
L_hcb = cell(1,nSubs);           

T_hcb = cell(1,nSubs);

Ms_rearranged = cell(1,nSubs);
Ks_rearranged = cell(1,nSubs);

VM = cell(1,nSubs);
CM = cell(1,nSubs);

for iSub = 1:nSubs
    if isempty(x)
        u0 = zeros(PrimalSub.Substructures(iSub).Mesh.nDOFs,1);
    else
        u0 = x{iSub};
    end
    
    %if the substructures have different number of modes
    if iscell(nM)
        iSubnM = nM{iSub};
    else
        iSubnM = nM;
    end
    
    Ms = PrimalSub.Substructures(iSub).mass_matrix();
    [Ks,~] = PrimalSub.Substructures(iSub).tangent_stiffness_and_force(u0);
    
    % Internal Mass and Stiffness matrix (constrained to
    % internal DOFs)
    Ms_II = Ms(PrimalSub.InternalFreeDOF{iSub},PrimalSub.InternalFreeDOF{iSub});
    Ks_II = Ks(PrimalSub.InternalFreeDOF{iSub},PrimalSub.InternalFreeDOF{iSub});
    
    
    [Phi,D] = eigs(Ks_II,Ms_II,iSubnM,'SM');
    [f0,ind] = sort(sqrt(diag(D))/2/pi);
    VMs = Phi(:,ind);
    
    % check on imaginary modes
    if ~isreal(f0)
        warning(' Complex eigenfrequencies')
        disp(' ')
    end
    %VMs = ROBnormalization(VMs,'mass',0,Ms_II);
    VM{iSub} = VMs;
    
    %static constraint modes
    Ks_IB = Ks(PrimalSub.InternalFreeDOF{iSub},PrimalSub.InterfaceDOF{iSub});  %self.BoundaryDOF{i}
    CMs = -Ks_II\Ks_IB;
    Us_internal = PrimalSub.InternalFreeDOF{iSub};
    Us_interface = PrimalSub.InterfaceDOF{iSub};
    Ms_rearranged{iSub} = [Ms(Us_interface,Us_internal) Ms(Us_interface,Us_interface);Ms_II Ms(Us_internal,Us_interface)];
    Ks_rearranged{iSub} = [Ks(Us_interface,Us_internal) Ks(Us_interface,Us_interface); Ks_II Ks_IB];
    
    %CMs = ROBnormalization(CMs,'displacement');
    CM{iSub} = CMs;
    
    B = length(PrimalSub.InterfaceDOF{iSub});
    
    %Reduction Basis
    T_hcb{iSub} = [eye(B)   zeros(B,size(VMs,2)) ;...
        CMs      VMs ];
    
end


nDOFPerNode = PrimalSub.Substructures(1).Mesh.nDOFPerNode;
nGlobalInterfaceDOFs = PrimalSub.nInt*nDOFPerNode;

%The currentGlobalDOF is used to compute the internal reduced
%DOFs and updated during the iteration trough the substructures
currentGlobalDOF = nGlobalInterfaceDOFs+1;

for iSub = 1:nSubs
    
    % Extract non-zero Interface nodes for Substructure $\texttt{i}$ from the
    % $\texttt{Interfaces}$ connectivity matrix.
    iSubInterface = PrimalSub.Interfaces(:,iSub);
    iSubNonZeroInterfaceNodes = logical(iSubInterface);
    
    % All Global Interface nodes
    AllGlobalInterfaceNodes = (1:PrimalSub.nInt)';
    % Obtain the unique nodes belonging to the interface.
    iSubGlobalInterfaceNodes = AllGlobalInterfaceNodes(iSubNonZeroInterfaceNodes);
    nInterfaceNodes = length(iSubGlobalInterfaceNodes);
    
    % The DOFs associated to a given Node $\texttt{n}$ in the local numbering
    % of a substructure are given by $\texttt{(n-1)*nDOFPerNode+1 : n*nDOFPerNode}$.
    iSubInterfaceDOFs = zeros(nInterfaceNodes*nDOFPerNode,1);
    for i = 1:nInterfaceNodes
        for j = 1:nDOFPerNode
            iSubInterfaceDOFs((i-1)*nDOFPerNode+j) = (iSubGlobalInterfaceNodes(i)-1)*nDOFPerNode+j;
        end
    end
    % store global Interface DOFs of the current Sub in the
    % global DOF vector
    Us_hcb{iSub} = iSubInterfaceDOFs;
    
    % Next we compute the reduced global Internal DOFs:
    iSubNreducedDOFs = size(T_hcb{iSub},2);
    NumberOfReducedInternalDOF = iSubNreducedDOFs - length(Us_hcb{iSub});
    
    % compute Global Reduced Internal DOFs in dependence of the
    % current Global DOFs which is precomputed before the loop
    % for the first iteration and afterwards updated
    iSubReducedInternalDOFs = (currentGlobalDOF : currentGlobalDOF+NumberOfReducedInternalDOF-1)';
    
    
    Us_hcb{iSub} = [Us_hcb{iSub}; iSubReducedInternalDOFs ];
    
    currentGlobalDOF = max(Us_hcb{iSub})+1;

end

nDOFg = max(Us_hcb{nSubs});

for iSub = 1: nSubs
    
    us = Us_hcb{iSub}; %reduced DOF vector of substructure i
    n_s = length(us);
    
    Ls = sparse(1:n_s, us, true(n_s,1), n_s, nDOFg );
    L_hcb{iSub} = Ls;
    
end

M_hcb = [];
K_hcb = [];

for iSub=1:nSubs
    %Ms_hcb = T_hcb{iSub}'*Ms_rearranged{iSub}*T_hcb{iSub};
    %Ks_hcb = T_hcb{iSub}'*Ks_rearranged{iSub}*T_hcb{iSub};
    Ms = PrimalSub.Substructures(iSub).mass_matrix();
    
    if isempty(x)
        u0 = zeros(PrimalSub.Substructures(iSub).Mesh.nDOFs,1);
    else
        u0 = x{iSub};
    end
    
    [Ks,~] = PrimalSub.Substructures(iSub).tangent_stiffness_and_force(u0);
    
    Mcs = PrimalSub.Substructures(iSub).constrain_matrix(Ms);
    Kcs = PrimalSub.Substructures(iSub).constrain_matrix(Ks);
    
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

