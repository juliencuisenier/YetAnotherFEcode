function [DS, names] = defect_sensitivities_substurcturing(PrimalSub, Phi, U, formulation, USEJULIA)


if nargin < 5
    formulation = 'N1';
    USEJULIA = 0;
    fprintf(' Default formulation for DS is: N1')
elseif nargin < 6
    USEJULIA = 0;
end

n  = size(Phi, 1);
nm = size(Phi, 2);
nd = size(U,   2);


K0 = PrimalSub.DATA.Kc;

DS = zeros(n, nm*nd);
names = zeros(nm*nd, 2);
cc = 1;
for kk = 1 : nd
    Uk = U(:,kk);  
    if USEJULIA == 1
        Ls = PrimalSub.L{1};
        Uks = Ls*Uk;
        dKs_dxi_k = stiffness_matrix_sensitivity(PrimalSub.Substructures(1), PrimalSub.Elements{1}, Uks, formulation);
        dK_dxi_k = Ls'*dKs_dxi_k*Ls;
        
        for iSub = 2:PrimalSub.nSubs
            Ls = PrimalSub.L{iSub};
            Uks = Ls*Uk;
            dKs_dxi_k = stiffness_matrix_sensitivity(PrimalSub.Substructures(iSub), PrimalSub.Elements{iSub}, Uks, formulation);
            dK_dxi_k = dK_dxi_k + Ls'*dKs_dxi_k*Ls;
        end
        
    else
        
        % Initializing dK_dxi_k with the 1rst substructure
        fprintf(' dKdxi(1), assembling %d elements ...', size(PrimalSub.Elements{1},1))
        t0 = tic;
        Ls = PrimalSub.L{1};
        Uks = Ls*Uk;
        dKs_dxi_k = PrimalSub.Substructures(1).matrix('stiffness_defect_derivative',Uks, formulation);
        dK_dxi_k = Ls'*dKs_dxi_k*Ls;
        fprintf(' %.2f s\n',toc(t0))
        
 
        % Adding the contribution of others substructures
        for iSub = 2:PrimalSub.nSubs
            fprintf(' dKdxi(%d), assembling %d elements ...',iSub, size(PrimalSub.Elements{iSub},1))
            t0 = tic;
            Ls = PrimalSub.L{iSub};
            Uks = Ls*Uk;
            dKs_dxi_k = PrimalSub.Substructures(iSub).matrix('stiffness_defect_derivative',Uks, formulation);
            dK_dxi_k = dK_dxi_k + Ls'*dKs_dxi_k*Ls;
            fprintf(' %.2f s\n',toc(t0))
        end
    end
    dK_dxi_k = PrimalSub.constrain_matrix( dK_dxi_k );
    for ii = 1 : nm
        Phi_i = PrimalSub.constrain_vector( Phi(:, ii) );
        Xi = - K0 \ (dK_dxi_k * Phi_i);
        DS(:, cc) =  PrimalSub.unconstrain_vector( Xi ) / max(abs(Xi));
        names(cc, :) = [ii kk]';
        cc = cc+1;
    end
end
disp(' ')
end

