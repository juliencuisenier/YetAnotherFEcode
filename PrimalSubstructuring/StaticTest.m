% Test of the PrimalSubstructuring class on a beam
clear
close all; 
clc

    
%% MODELs (material, mesh, assemblies)                              

 %elementType = 'HEX20';
elementType = 'HEX8';

% PREPARE MODEL
% MATERIAL ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 

myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

% Element

myElementConstructor = @()Hex8Element(myMaterial);


% SUBMESHING:______________________________________________________________

TotalMesh = load_gmsh2("RotorSketch2.msh", -1); %reading the .msh file

Submeshes = Submeshing(TotalMesh);


%Mesh 1 (Left blade)_______________________________________________________

myMesh1 = Mesh(Submeshes{1,1});
myMesh1.create_elements_table(Submeshes{1,2}, myElementConstructor);

%Mesh 2 (Top blade)________________________________________________________

myMesh2 = Mesh(Submeshes{2,1});
myMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);

%Mesh 3 (Right blade)______________________________________________________

myMesh3 = Mesh(Submeshes{3,1});
myMesh3.create_elements_table(Submeshes{3,2}, myElementConstructor);

%Mesh 4 (Bottom blade)_____________________________________________________

myMesh4 = Mesh(Submeshes{4,1});
myMesh4.create_elements_table(Submeshes{4,2}, myElementConstructor);

%Mesh 5 (Heart)____________________________________________________________

myMesh5 = Mesh(Submeshes{5,1});
myMesh5.create_elements_table(Submeshes{5,2}, myElementConstructor);
myMesh5.set_essential_boundary_condition(Submeshes{5,3}{1},1:3,0)
myMesh5.set_essential_boundary_condition(Submeshes{5,3}{4},1:3,0)



% ASSEMBLY ________________________________________________________________
% ReducedAssembly is a subclass of Assembly. It adds methods for the
% computation of reduced internal forces, tangent stiffness matrix,
% stiffness tensors (for nominal and defected structures) and to project
% matrices and vectors in general using the basis V (new property of the
% class). Additionally, the new property "U" can be used to add a basis of
% defects to the structure.
% [Notice that ReducedAssembly inherits all properties and methods of the
% standard Assembly class]
Assembly1 = Assembly(myMesh1);
Assembly2 = Assembly(myMesh2);
Assembly3 = Assembly(myMesh3);
Assembly4 = Assembly(myMesh4);
Assembly5 = Assembly(myMesh5);

%% SUBSTRUCTURING                                                   

% PrimalSubstructuring is a class with the basics of the substructuring and
% uses the primal assembly model of resolution 

globalIndices = get_globalIndices(Submeshes);
PrimalSub = PrimalSubstructuring([Assembly1 Assembly2...
    Assembly3 Assembly4 Assembly5],globalIndices,[]);


%%REFERENCE MODEL__________________________________________________________
[nodes_ref,elements_ref,nset_ref] = extract_gmsh(TotalMesh);

Mesh_ref = Mesh(nodes_ref);
Mesh_ref.create_elements_table(elements_ref, myElementConstructor)

%Boundaries conditions

%At front of the heart
x_front = find(nodes_ref(:,1)==1);
y_bc = find(abs(nodes_ref(:,2))<=0.5);
z_bc = find(abs(nodes_ref(:,3))<=0.5);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);

Mesh_ref.set_essential_boundary_condition(indices_front,1:3,0)

%At behind of the heart

x_behind = find(nodes_ref(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

Mesh_ref.set_essential_boundary_condition(indices_behind,1:3,0)


%Assembly

Assembly_ref = Assembly(Mesh_ref);

M_ref = Assembly_ref.mass_matrix();
u0 = zeros( Mesh_ref.nDOFs, 1);
[K_ref,~] = Assembly_ref.tangent_stiffness_and_force(u0);

Mc_ref = Assembly_ref.constrain_matrix(M_ref);
Kc_ref = Assembly_ref.constrain_matrix(K_ref);

figure
hold on
PlotMesh(nodes_ref, elements_ref, 0);
legend('Mesh ref')

%% STATIC RESOLUTION

%With substructuring

%Sub1
f1=zeros(myMesh1.nDOFs,1);
f1(98) = -1e-6; %Force applied on the 2nd DOF of the 33rd node

%Sub2
f2=zeros(myMesh2.nDOFs,1);
f2(96) = 1e-6; %Force applied on the 3th DOF of the 32nd node

%Sub3
f3=zeros(myMesh3.nDOFs,1);
f3(65) = 1e-6; %Force applied on the 2nd DOF of the 22nd node

%Sub4
f4=zeros(myMesh4.nDOFs,1);
f4(36) = -1e-6; %Force applied on the 3rd DOF of the 12th node

%Sub5
f5=zeros(myMesh5.nDOFs,1);

Fext = {f1,f2,f3,f4,f5};

u = PrimalSub.static_resolution([],Fext);

%With global reference model

Fext_ref = zeros(Mesh_ref.nDOFs,1);

Fext_ref(254) = -1e-6; %Force of Sub1
Fext_ref(168) = 1e-6; %Force of Sub2
Fext_ref(392) = 1e-6; %Force of Sub3
Fext_ref(279) = -1e-6; %Force of Sub4

nNodes = size(nodes_ref,1);
u0 = zeros( Mesh_ref.nDOFs, 1);
[K,~] = Assembly_ref.tangent_stiffness_and_force(u0);

v = K\Fext_ref;

