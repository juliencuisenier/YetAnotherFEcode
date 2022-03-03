% Test of the PrimalSubstructuring class on a beam
clear
close all; 
clc

    
%% MODELs (material, mesh, assemblies)                              

% elementType = 'HEX20';
elementType = 'TET10';

% PREPARE MODEL
% MATERIAL ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

% Element

myElementConstructor = @()Tet10Element(myMaterial);

% MESH:__________________________________________________________________

TotalMesh = load_gmsh2("GeomBeam.msh", 11); %reading the .msh file

[Substructures] = Submeshing(TotalMesh,2);

Substructures = reindexing_elements_from_global(Substructures);

Substructures = find_nsets(Substructures);

% ASSEMBLY ________________________________________________________________




% %% SUBSTRUCTURING                                                   
% 
% % PrimalSubstructuring is a class with the basics of the substructuring and
% % uses the primal assembly model of resolution 
% PrimalSub = PrimalSubstructuring([Assembly1 Assembly2]);
% 
% % INTERFACE _______________________________________________________________
% % Specify how the substructures are connected.
% 
% NrSub1 = 1;                                     % substructure #1
% NrSub2 = 2;                                     % substructure #2
% Inodes1 = nset1{4};                             % interface nodes of ss#1
% Inodes2 = nset2{1};                             % interface nodes of ss#2
% Interface_node1 = find_node(1.5,0,0,nodes1);	% node of ss#1 matching ...
% Interface_node2 = find_node(0,0,0,nodes2);      % ... a node of ss#2
% 
% %only used one edge or boundary point for the interface connection as there
% % exist no symetrise on the edge point (connection is only in one direction
% % possible). When errors occur then the an additional edge point should be
% % given as input
% create_Interface(PrimalSub, NrSub1, NrSub2, Inodes1, Inodes2, ...
%     Interface_node1, Interface_node2);  % --> property "Interfaces" added
%                                         %     to DualSub
%                                         % DualSub.Interfaces is a table 
%                                         % with the nodes of the interface 
%                                         % for ss#1 (first column) and ss#2
%                                         % (second column)
% 
% % GLOBAL COORDINATES ______________________________________________________
% % Translate and/or rotates the nodes of the substructures to make them
% % compatible. The shift of the coordinates is needed only for plotting
% % reasons. The fields "PrimalSub.Substructures(s).Mesh.Nodes" are updated
% % accordingly.
% transform_substructures(PrimalSub);
% nodes2 = Assembly2.Mesh.nodes;      % update nodes of ss#2
% 
% 
% % Test of the static resolution ___________________________________________
% 
% % Some arbitrary external forces, fextN corresponding to the external force
% % applied on SubN
% fext1 = zeros(417,1);
% fext1(1:10)= 0.005*ones(10,1);
% 
% fext2 = zeros(573,1);
% 
% Fext = {fext1,fext2};
% 
% % Test of the static_resolution method by comparing the global resolution
% % with the substuctured one
% 
% % u = PrimalSub.global_static_resolution([],Fext);
% % 
% % v = PrimalSub.local_static_resolution([],Fext);
% 
% %% VIBRATION MODES                                                  
% PrimalSub.localization_matrix()
% PrimalSub.compute_Dirichlet_and_global_DOFs()
% % GLOBAL Mass and Stiffness matrices ______________________________________
% % Example: M = \sum_s(Ls'*Ms*Ls)
% [M,K] = PrimalSub.global_mass_stiffness([]);  % free
% Mc = PrimalSub.constrain_matrix(M);           % constrained
% Kc = PrimalSub.constrain_matrix(K);           % constrained
% 
% % EIGENMODES ______________________________________________________________
% n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
% 
% [V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
% [f0,ind] = sort(sqrt(diag(om))/2/pi);
% V0 = V0(:,ind);
% f0_full = f0;
% 
% 
%  
% V0  = PrimalSub.unconstrain_vector(V0);
% %V0  = ROBnormalization( V0, 'displacement' );
% V0s = deformation2substructs(PrimalSub,V0);   % Localize V0 to substructure s
% 
% % PLOT of the two substructures ___________________________________________
% mod = 2;
% figure
% PlotMesh(nodes1, elements1, 0);
% PlotMesh(nodes2, elements2, 0);
% nodalDef1 = reshape(V0s{1}(:,mod),3,[]).';
% nodalDef2 = reshape(V0s{2}(:,mod),3,[]).';
% PlotFieldonDeformedMesh(nodes1, elements1, nodalDef1, 'factor', 0.1)
% PlotFieldonDeformedMesh(nodes2, elements2, nodalDef2, 'factor', 0.1)
% colormap jet
% title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ...
%     ' Hz with substructuring'])