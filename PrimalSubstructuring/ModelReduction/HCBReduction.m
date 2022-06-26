classdef HCBReduction < PrimalSubstructuring
    %HCBREDUCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        nDOFglobal_hcb
        ms %array containing the reduced internal DOFs
        Us_hcb
        T_hcb %cell of the projection matrices of each substructure
        L_hcb %localization matrices
    end
    
    methods
        function self = HCBReduction(PrimalSub,freq)
            self.nSubs = PrimalSub.nSubs;
            self.Interface = PrimalSub.Interfaces
        end
        
        function T_hcb = compute_T(self,PrimalSub,freq)
            
            nSubs = self.nSubs;
            self.ms = zeros(1,nSubs);
            self.Us_hcb = cell(1,nSubs);
            self.T_hcb = cell(1,nSubs);
            
  
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
                
                self.T_hcb{iSub} = [IVM CM]; 
                self.ms(iSub) = m;
                
                self.Us_hcb{iSub} = zeros(1,m+nInts);    
            end
        end
        
        function nDOFglobal_hcb = compute_global_DOF_hcb(self,
        
    end
end

