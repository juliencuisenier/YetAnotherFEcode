% Test of the PrimalSubstructuring class on a beam
clear
close all; 
clc

    
%% MODELs (material, mesh, assemblies)                              

% elementType = 'HEX20';
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

% SUBMESHING:__________________________________________________________________

TotalMesh = load_gmsh2("GeomBeam.msh", 5); %reading the .msh file

Submeshes = Submeshing(TotalMesh,2);

Submeshes = reindexing_elements_from_global(Submeshes);

Submeshes = find_nsets(Submeshes);

%Mesh 1___________________________________________________________________

myMesh1 = Mesh(Submeshes{1,1});
myMesh1.create_elements_table(Submeshes{1,2}, myElementConstructor);
myMesh1.set_essential_boundary_condition(Submeshes{1,4}{5},1:3,0)

figure
hold on
PlotMesh(Submeshes{1,1}, Submeshes{1,2}, 0);
legend('Mesh SS#1')

%Mesh 2____________________________________________________________________

myMesh2 = Mesh(Submeshes{2,1});
myMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);
myMesh2.set_essential_boundary_condition(Submeshes{2,4}{2},1:3,0)

figure
hold on
PlotMesh(Submeshes{2,1}, Submeshes{2,2}, 0);
legend('Mesh SS#2')

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

%% SUBSTRUCTURING                                                   

% PrimalSubstructuring is a class with the basics of the substructuring and
% uses the primal assembly model of resolution 
PrimalSub = PrimalSubstructuring([Assembly1 Assembly2]);

%Interface

create_Interface(PrimalSub,Submeshes,1, 2)


%%REFERENCE MODEL__________________________________________________________
[nodes_ref,elements_ref,nset_ref] = extract_gmsh(TotalMesh);

Mesh_ref = Mesh(nodes_ref);
Mesh_ref.create_elements_table(elements_ref, myElementConstructor)
Mesh_ref.set_essential_boundary_condition(nset_ref{2},1:3,0)
Mesh_ref.set_essential_boundary_condition(nset_ref{5},1:3,0)

Assembly_ref = Assembly(Mesh_ref);

figure
hold on
PlotMesh(nodes_ref, elements_ref, 0);
legend('Mesh ref')

%%STATIC RESOLUTIONS_______________________________________________________

%With substructuring
fext1 = zeros(3699,1);
fext1(250:260)=1;

fext2 = zeros(3699,1);

Fext = {fext1,fext2};

u = PrimalSub.static_resolution([],Fext);

%With global reference model

Fext_ref = zeros(7203,1);
Fext_ref(Submeshes{1,3}(250:260)) = 1;

nNodes = size(nodes_ref,1);
u0 = zeros( Mesh_ref.nDOFs, 1);
[K,~] = Assembly_ref.tangent_stiffness_and_force(u0);

v = K\Fext_ref;
