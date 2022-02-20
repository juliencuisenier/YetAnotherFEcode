classdef PrimalSubstructuring < handle
    properties
        Substructures    % Array of Assembly class objects: these are the substructures
                         % The first substructure in the array is the
                         % master substructure
        nSubs            % Total number of substructures
        
        Interfaces = []  % matrix of interface connectivity (nInt X xSubs) (Nodes)
        nInt             % Total number of Interface connected Nodes
        InterfaceDOF = {}% cell array containing the Interface DOFs for each substructure
        InternalFreeDOF = {} % cell array containing the Internal free DOFs for each substructure
        
        nDOFglobal       %global number of DOFs
        globalFreeDOFs   %global free DOFs
        B                %sparse global boundary matrix
        DirichletDOFs = []%two column vector first global DOF second dirichlet value
        
        Us = {}          %cell array with global DOF indices vectors of the substructures
        L = {}           %Localization matrix
        
        DATA             %Miscellaneous user-defined data which can be stored
                         %e.g. DATA.K, DATA.M, DATA.C etc.
    end
    
    methods
        function self = PrimalSubstructuring(Substructures)
           self.Substructures = Substructures;
           self.nSubs = self.nSubs();
        end 
        
        function nSubs = get.nSubs(self)
            nSubs = length(self.Substructures);
        end        
  
        function nInt = get.nInt(self)
            nInt = size(self.Interfaces,1);
        end 
        
        function localization_matrix(self)

            %for explenation of the procedure see Method 2 in example_localization.mlx
            self.compute_InterfaceDOF();
            self.compute_global_DOF();
            nDOFg = self.nDOFglobal; %global DOF vector
            
            for i = 1: self.nSubs
                
                us = self.Us{i}; %DOF vector of substructure i
                n_s = length(us);

                Ls = sparse(1:n_s, us, true(n_s,1), n_s, nDOFg );
                self.L{i} = Ls;

            end
        end
        
        function compute_InterfaceDOF(self)
                %This functions computes the Interface DOFs for all
                %substructures. A condition is that all substructures have the
                %same nDOFs per Node.
                nDOFPerNode = self.Substructures(1).Mesh.nDOFPerNode;   

                for iSub = 1:self.nSubs            
                % Extract non-zero Interface nodes for Substructure $\texttt{i}$ from the 
                % $\texttt{Interfaces}$ connectivity matrix. 
                    iSubInterface = self.Interfaces(:,iSub);
                    iSubNonZeroInterfaceNodes = logical(iSubInterface);    

                % Obtain the unique nodes belonging to the interface. 
                    iSubInterfaceNodes = self.Interfaces(iSubNonZeroInterfaceNodes,iSub);
                    nInterfaceNodes = length(iSubInterfaceNodes);    
                % The DOFs associated to a given Node $\texttt{n}$ in the local numbering 
                % of a substructure are given by $\texttt{(n-1)*nDOFPerNode+1 : n*nDOFPerNode}$. 
                    iSubInterfaceDOFs = zeros(nInterfaceNodes*nDOFPerNode,1);
                    for i = 1:nInterfaceNodes
                        for j = 1:nDOFPerNode
                            iSubInterfaceDOFs((i-1)*nDOFPerNode+j) = (iSubInterfaceNodes(i)-1)*nDOFPerNode+j;

                        end
                    end
                    self.InterfaceDOF{iSub} = iSubInterfaceDOFs;
                end

        end
        
        function compute_internal_freeDOF(self)
            self.InternalFreeDOF = cell(1,self.nSubs);
            for i = 1:self.nSubs
                self.InternalFreeDOF{i} = ...
                    setdiff((1 : self.Substructures(i).Mesh.nDOFs)', self.InterfaceDOF{i});
                
                if ~isempty(self.Substructures(i).Mesh.BC.Dirichlet)
                    self.InternalFreeDOF{i} = ...
                        setdiff(self.InternalFreeDOF{i}, self.Substructures(i).Mesh.BC.Dirichlet(:,1));
                end                  
            end
            
            
        end
        
        function compute_global_DOF(self)
            %
            % This method computes the global DOFs for each
            % substructure and stores it in the cell array Us{iSub}
            %In this function it is assumed that all substructures have
            %the same nDOFPerNode...
            nDOFPerNode = self.Substructures(1).Mesh.nDOFPerNode;
            nGlobalInterfaceDOFs = self.nInt*nDOFPerNode;
            
            %The currentGlobalDOF is used to compute the internal reduced
            %DOFs and updated during the iteration trough the substructures
            currentGlobalDOF = nGlobalInterfaceDOFs+1;
            
            for iSub = 1:self.nSubs
                
                self.Us{iSub} = zeros(self.Substructures(iSub).Mesh.nDOFs,1);
            % Extract non-zero Interface nodes for Substructure $\texttt{i}$ from the 
            % $\texttt{Interfaces}$ connectivity matrix. 
                
                iSubInterface = self.Interfaces(:,iSub);
                iSubNonZeroInterfaceNodes = logical(iSubInterface);    
                
            % All Global Interface nodes
                AllGlobalInterfaceNodes = (1:self.nInt)';
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
                self.Us{iSub}(self.InterfaceDOF{iSub}) = iSubInterfaceDOFs;
            
                % Next we compute the global Internal DOFs:
                localNonInterfaceDOFs = setdiff((1 : self.Substructures(iSub).Mesh.nDOFs)', self.InterfaceDOF{iSub});
                
                iSubHighestGlobalDOF = ...
                    length(localNonInterfaceDOFs) + currentGlobalDOF -1;
                
                self.Us{iSub}(localNonInterfaceDOFs) = ...
                    (currentGlobalDOF : iSubHighestGlobalDOF)';
                
                currentGlobalDOF = iSubHighestGlobalDOF + 1;
            end
            
        end
        
        function nDOFglobal = get.nDOFglobal(self)
            nDOFglobal = max(self.Us{self.nSubs});
        end 
        
        
        function compute_Dirichlet_and_global_DOFs(self)
        % rearrange all dirichlet BC from the Assemblies in global
        % Dirichlet BC and compute the global_FreeDOFs and the boolean
        % matrix B which is used to uncostrain vectors (same principle as
        % for the normal assembly but now on a global level)
            for iSub = 1:self.nSubs
             
                 iUs = self.Us{iSub};
                 
                 %check each substructure if there are Dirichlet BCs and
                 %add them to the global DirichletDOFs
                 if ~isempty(self.Substructures(iSub).Mesh.EBC.constrainedDOFs)
                     DirichletDOFsLocal = self.Substructures(iSub).Mesh.EBC.constrainedDOFs;
                     DirichletLocal = [iUs(DirichletDOFsLocal(:,1)) DirichletDOFsLocal(:,2)];
                     self.DirichletDOFs = [self.DirichletDOFs; DirichletLocal];
                 end
            end
            
            self.globalFreeDOFs = setdiff(1 : self.nDOFglobal, self.DirichletDOFs);
            
            % update the boolean matrix
            nf = length(self.globalFreeDOFs);
            self.B = sparse(self.globalFreeDOFs,1:nf,true,self.nDOFglobal,nf);
            
        end
        
        function [M,K] = global_mass_stiffness(self,x)
            %Global_Mass_Stiffness matrix
            %Computes the global mass and stiffness matrix of all
            %substructures
            %Input:
            %x --> system data (e.g. initial deformation)
            M = [];
            K = [];
            for iSub = 1:self.nSubs
                 
                if isempty(x)
                    u0 = zeros(length(self.Us{iSub}),1);
                else
                    u0 = x{iSub};
                end

                 Ls = self.L{iSub};
                 

                 [Ks,~]=self.Substructures(iSub).tangent_stiffness_and_force(u0);

                 Ms = self.Substructures(iSub).mass_matrix();
                
                 %With the localization matrix we can compute the global
                 %matrices
                 Ms = Ls'*Ms*Ls;
                 Ks = Ls'*Ks*Ls;
                 
                 if isempty(M)
                    M = Ms;
                    K = Ks;
                 else
                    M = M + Ms;
                    K = K + Ks;
                 end
            end
        end
        
        
        function Mc = constrain_matrix(self,Matrix)
            
            Mc = Matrix(self.globalFreeDOFs,self.globalFreeDOFs);
            
        end
        
        function vc = constrain_vector(self,v)
            
            vc = v(self.globalFreeDOFs,:);
            
        end
        
        function v = unconstrain_vector(self,vc)
            v = self.B * vc;
            v(self.DirichletDOFs(:,1),:) = repmat(self.DirichletDOFs(:,2),1,size(vc,2));
        end
        
        function fg = ps_force(self,Fext) %where Fext is a cell of the
                                          %external forces applied on each
           fg = [];                       %substuctures
           
           for iSub = 1:self.nSubs 
               Ls = self.L{iSub};
               
               if isempty(fg)
                   fg = Ls'*Fext{iSub};
               else
                   fg = fg + Ls'*Fext{iSub};
               end
           end
            
        end
        
        function u = global_static_resolution(self,x,Fext)
            self.localization_matrix()
            [M,K] = self.global_mass_stiffness(x);
            fg = self.ps_force(Fext); 
            u = (M+K)\fg;
        end
        
        function u = local_static_resolution(self,x,Fext)
            self.localization_matrix()
            
            u = [];
            
            for iSub=1:self.nSubs
                if isempty(x)
                    u0 = zeros(length(self.Us{iSub}),1);
                else
                    u0 = x{iSub};
                end
                    
                [Ks,~]=self.Substructures(iSub).tangent_stiffness_and_force(u0);

                 Ms = self.Substructures(iSub).mass_matrix();
                 
                 us{iSub} = (Ms+Ks)\Fext{iSub};
                 
                 Ls = self.L{iSub};
                 
                 if isempty(u)
                     u = Ls'*us{iSub};
                 else
                     u = u + Ls'*us{iSub};
                 end
            end    
        end
        
        
    end
    

end