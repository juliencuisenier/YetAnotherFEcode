% Test of the PrimalSubstructuring class on a beam
close all; 
clc

    
%% MODELs (material, mesh, assemblies)                              

% elementType = 'HEX20';
elementType = 'TET10';

% PREPARE MODEL
% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
myMaterial = KirchoffMaterial(E,rho,nu);

% Element
switch elementType
    case 'HEX20'
        myelement = Hex20DefectsElement(myMaterial);
    case 'TET10'
        myelement = Tet10DefectsElement(myMaterial);
end

% MESH_1:__________________________________________________________________
l1 = 0.2;
w1 = 0.1;
t1 = 0.008; % 0.015;
nx1 = 5;
ny1 = 1;
nz1 = 1;
[nodes1, elements1, nset1] = ...
    mesh_3Dparallelepiped(elementType, l1, w1, t1, nx1, ny1, nz1);
nDOFperNode = 3;
myMesh1 = Mesh(nodes1, elements1);
myMesh1.set_nDOFPerNode(nDOFperNode);
myMesh1.ElementsTable = create_elements_table(elements1, myelement);

% MESH > boundary conditions
myMesh1.BC.set_dirichlet_dofs(nset1{1}, 1:3, 0)

figure
hold on
PlotMesh(nodes1, elements1, 0);
nodeplot(nodes1, nset1{1}, 'r') % BCs
nodeplot(nodes1, nset1{4}, 'c') % Interface
legend('Mesh SS#1', 'BCs', 'Interface')

% MESH_2:__________________________________________________________________
l2 = 0.3;
w2 = 0.1;
t2 =  0.008; % 0.015;
nx2 = 7;
ny2 = 1;
nz2 = 1;
[nodes2, elements2, nset2]= ...
    mesh_3Dparallelepiped(elementType, l2, w2, t2, nx2, ny2, nz2);

nDOFperNode = 3;
myMesh2 = Mesh(nodes2, elements2);
myMesh2.set_nDOFPerNode(nDOFperNode);
myMesh2.ElementsTable = create_elements_table(elements2, myelement);

% MESH > boundary conditions
myMesh2.BC.set_dirichlet_dofs(nset2{4}, 1:3, 0)

figure
hold on
PlotMesh(nodes2 ,elements2, 0);
nodeplot(nodes2, nset2{4}, 'r') % BCs
nodeplot(nodes2, nset2{1}, 'c') % Interface
legend('Mesh SS#2', 'BCs', 'Interface')


% ASSEMBLY ________________________________________________________________
% ReducedAssembly is a subclass of Assembly. It adds methods for the
% computation of reduced internal forces, tangent stiffness matrix,
% stiffness tensors (for nominal and defected structures) and to project
% matrices and vectors in general using the basis V (new property of the
% class). Additionally, the new property "U" can be used to add a basis of
% defects to the structure.
% [Notice that ReducedAssembly inherits all properties and methods of the
% standard Assembly class]
Assembly1 = ReducedAssembly(myMesh1);
Assembly2 = ReducedAssembly(myMesh2);


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
% reasons. The fields "DualSub.Substructures(s).Mesh.Nodes" are updated
% accordingly.
transform_substructures(PrimalSub);
nodes2 = Assembly2.Mesh.Nodes;      % update nodes of ss#2


% Test of the static resolution ___________________________________________

% Some arbitrary external forces, fextN corresponding to the external force
% applied on SubN
fext1 = zeros(417,1);
fext1(1:50) = 50*ones(50,1);

fext2 = zeros(573,1);

Fext = {fext1,fext2};

% Test of the static_resolution method
u = PrimalSub.static_resolution([],Fext);
