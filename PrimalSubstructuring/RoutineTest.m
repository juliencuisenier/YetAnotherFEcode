% Test of the PrimalSubstructuring class on a beam
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

% MESH_1:__________________________________________________________________
l1 = 0.2;
w1 = 0.1;
t1 = 0.008; % 0.015;
nx1 = 5;
ny1 = 1;
nz1 = 1;
[nodes1, elements1, nset1] = ...
    mesh_3Dparallelepiped(elementType, l1, w1, t1, nx1, ny1, nz1);

myMesh1 = Mesh(nodes1);
myMesh1.create_elements_table(elements1, myElementConstructor);

% MESH > boundary conditions
%myMesh1.BC.set_dirichlet_dofs(nset1{1}, 1:3, 0)
myMesh1.set_essential_boundary_condition(nset1{1},1:3,0)

figure
hold on
PlotMesh(nodes1, elements1, 0);
legend('Mesh SS#1')

% MESH_2:__________________________________________________________________
l2 = 0.3;
w2 = 0.1;
t2 =  0.008; % 0.015;
nx2 = 7;
ny2 = 1;
nz2 = 1;
[nodes2, elements2, nset2]= ...
    mesh_3Dparallelepiped(elementType, l2, w2, t2, nx2, ny2, nz2);


myMesh2 = Mesh(nodes2);
myMesh2.create_elements_table(elements2, myElementConstructor);

% MESH > boundary conditions
%myMesh2.BC.set_dirichlet_dofs(nset2{4}, 1:3, 0)
myMesh2.set_essential_boundary_condition(nset2{4},1:3,0)

figure
hold on
PlotMesh(nodes2 ,elements2, 0);
legend('Mesh SS#2')


% ASSEMBLY ________________________________________________________________

Assembly1 = Assembly(myMesh1);
Assembly2 = Assembly(myMesh2);


%% SUBSTRUCTURING                                                   

% PrimalSubstructuring is a class with the basics of the substructuring and
% uses the primal assembly model of resolution 
PrimalSub = PrimalSubstructuring([Assembly1 Assembly2]);

% INTERFACE _______________________________________________________________
% Specify how the substructures are connected.

NrSub1 = 1;                                     % substructure #1
NrSub2 = 2;                                     % substructure #2
Inodes1 = nset1{4};                             % interface nodes of ss#1
Inodes2 = nset2{1};                             % interface nodes of ss#2
Interface_node1 = find_node(1.5,0,0,nodes1);	% node of ss#1 matching ...
Interface_node2 = find_node(0,0,0,nodes2);      % ... a node of ss#2

%only used one edge or boundary point for the interface connection as there
% exist no symetrise on the edge point (connection is only in one direction
% possible). When errors occur then the an additional edge point should be
% given as input
create_Interface(PrimalSub, NrSub1, NrSub2, Inodes1, Inodes2, ...
    Interface_node1, Interface_node2);  % --> property "Interfaces" added
                                        %     to DualSub
                                        % DualSub.Interfaces is a table 
                                        % with the nodes of the interface 
                                        % for ss#1 (first column) and ss#2
                                        % (second column)

% GLOBAL COORDINATES ______________________________________________________
% Translate and/or rotates the nodes of the substructures to make them
% compatible. The shift of the coordinates is needed only for plotting
% reasons. The fields "PrimalSub.Substructures(s).Mesh.Nodes" are updated
% accordingly.
transform_substructures(PrimalSub);
nodes2 = Assembly2.Mesh.nodes;      % update nodes of ss#2


% Test of the static resolution ___________________________________________

% Some arbitrary external forces, fextN corresponding to the external force
% applied on SubN
fext1 = zeros(417,1);
fext1(1:10)= 0.005*ones(10,1);

fext2 = zeros(573,1);

Fext = {fext1,fext2};

% Test of the static_resolution method by comparing the global resolution
% with the substuctured one

% u = PrimalSub.global_static_resolution([],Fext);
% 
% v = PrimalSub.local_static_resolution([],Fext);

%% VIBRATION MODES                                                  
PrimalSub.localization_matrix()
PrimalSub.compute_Dirichlet_and_global_DOFs()
% GLOBAL Mass and Stiffness matrices ______________________________________
% Example: M = \sum_s(Ls'*Ms*Ls)
[M,K] = PrimalSub.global_mass_stiffness([]);  % free
Mc = PrimalSub.constrain_matrix(M);           % constrained
Kc = PrimalSub.constrain_matrix(K);           % constrained

% EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

[V0,om] = eigs(Kc, Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);
f0_full = f0;


 
V0  = PrimalSub.unconstrain_vector(V0);
%V0  = ROBnormalization( V0, 'displacement' );
V0s = deformation2substructs(PrimalSub,V0);   % Localize V0 to substructure s

% PLOT of the two substructures ___________________________________________
mod = 2;
figure
PlotMesh(nodes1, elements1, 0);
PlotMesh(nodes2, elements2, 0);
nodalDef1 = reshape(V0s{1}(:,mod),3,[]).';
nodalDef2 = reshape(V0s{2}(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes1, elements1, nodalDef1, 'factor', 0.1)
PlotFieldonDeformedMesh(nodes2, elements2, nodalDef2, 'factor', 0.1)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ...
    ' Hz with substructuring'])
