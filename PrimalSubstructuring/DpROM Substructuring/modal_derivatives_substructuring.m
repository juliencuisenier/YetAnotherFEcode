function [MD, names] = modal_derivatives_substructuring(PrimalSub, Phi, USEJULIA)

if nargin < 4
    USEJULIA = 0;
end

n = size(Phi,1);
n_VMs = size(Phi,2);


K0 = PrimalSub.DATA.Kc;

MD = zeros(n, n_VMs*(n_VMs+1)/2);
names = zeros(n_VMs*(n_VMs+1)/2, 2);
kk = 1;
for jj = 1 : n_VMs
    
    Phi_j = Phi(:, jj);
    if USEJULIA == 1
        
        % Initializing dK_deta_j with the 1rst substructure
        Phis_j = PrimalSub.L{1}*Phi_j;
        dKs_deta_j = stiffness_matrix_derivative(PrimalSub.Substructures(1), PrimalSub.Elements{1}, Phis_j);
        dK_deta_j = PrimalSub.L{1}'*dKs_deta_j*PrimalSub.L{1};
        
        % Adding the contribution of others substructures
        for iSub = 2:PrimalSub.nSubs
            Phis_j = PrimalSub.L{iSub}*Phi_j;
            dKs_deta_j = stiffness_matrix_derivative(PrimalSub.Substructures(iSub), PrimalSub.Elements{iSub}, Phis_j);
            dK_deta_j = dK_deta_j + PrimalSub.L{iSub}'*dKs_deta_j*PrimalSub.L{iSub};
        end
        
    else
        % Initializing dK_deta_j with the 1rst substructure
        t0 = tic;
        fprintf(' dKsdq(1), asembling %d elements ...', size(elements,1))
        Phis_j = PrimalSub.L{1}*Phi_j;
        dKs_deta_j = PrimalSub.Substructures(1).matrix('stiffness_derivative',Phis_j);
        dK_deta_j = PrimalSub.L{1}'*dKs_deta_j*PrimalSub.L{1};
        fprintf(' %.2f s\n',toc(t0))
        
        % Adding the contribution of others substructures
        for iSub = 2:PrimalSub.nSubs
            t0 = tic;
            fprintf(' dKsdq(%d), asembling %d elements ...',iSub, size(PrimalSub.Elements{iSub},1))
            Phis_j = PrimalSub.L{iSub}*Phi_j;
            dKs_deta_j = PrimalSub.Substructures(iSub).matrix('stiffness_derivative',Phis_j);
            dK_deta_j = dK_deta_j + PrimalSub.L{iSub}'*dKs_deta_j*PrimalSub.L{iSub};
            fprintf(' %.2f s\n',toc(t0))
        end
    end
    dK_deta_j = PrimalSub.constrain_matrix( dK_deta_j );
    
    for ii = 1 : n_VMs
        if ii < jj
            continue
        end
        
        Phi_i = PrimalSub.constrain_vector( Phi(:, ii) );
        dPhi_i_deta_j = -K0\(dK_deta_j * Phi_i); 
        
        th =  dPhi_i_deta_j / max(abs(dPhi_i_deta_j));
        MD(:,kk) = PrimalSub.unconstrain_vector( th );
        names(kk, :) = [ii jj];
        kk = kk + 1;
    end
end
