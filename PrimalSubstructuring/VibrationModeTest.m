% Test of the PrimalSubstructuring class on a clambed-clambed beam
clear
close all; 
clc

    
%% MODELs (material, mesh, assemblies)                              

 elementType = 'HEX20';
%elementType = 'HEX8';

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

TotalMesh = load_gmsh2("geomBeam - Copie.msh", -1); %reading the .msh file

Submeshes = Submeshing(TotalMesh);


%Mesh 1___________________________________________________________________

myMesh1 = Mesh(Submeshes{1,1});
myMesh1.create_elements_table(Submeshes{1,2}, myElementConstructor);
myMesh1.set_essential_boundary_condition(Submeshes{1,3}{5},1:3,0)

%Mesh 2____________________________________________________________________

myMesh2 = Mesh(Submeshes{2,1});
myMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);
myMesh2.set_essential_boundary_condition(Submeshes{2,3}{2},1:3,0)


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
globalIndices = [Submeshes{1,4} Submeshes{2,4}];

PrimalSub = PrimalSubstructuring([Assembly1 Assembly2],globalIndices);

 
%% VIBRATION MODES                                                  

t1 = tic;
% EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

Mods = [2 4];

V0s = PrimalSub.vibration_mode(n_VMs,Mods);

t1f = toc(t1);

%REFERENCE MODEL 

t2 = tic;

mod = 2;

[nodes_ref,elements_ref,nset_ref] = extract_gmsh(TotalMesh);

Mesh_ref = Mesh(nodes_ref);
Mesh_ref.create_elements_table(elements_ref, myElementConstructor)
Mesh_ref.set_essential_boundary_condition(nset_ref{2},1:3,0)
Mesh_ref.set_essential_boundary_condition(nset_ref{5},1:3,0)

Assembly_ref = Assembly(Mesh_ref);

M_ref = Assembly_ref.mass_matrix();
u0 = zeros( Mesh_ref.nDOFs, 1);
[K_ref,~] = Assembly_ref.tangent_stiffness_and_force(u0);

Mc_ref = Assembly_ref.constrain_matrix(M_ref);
Kc_ref = Assembly_ref.constrain_matrix(K_ref);

[V0_ref,om_ref] = eigs(Kc_ref, Mc_ref, n_VMs, 'SM');
[f0_ref,ind_ref] = sort(sqrt(diag(om_ref))/2/pi);
V0_ref = V0_ref(:,ind_ref);

V0_ref = Assembly_ref.unconstrain_vector(V0_ref);

figure
hold on
nodalDef_ref = reshape(V0_ref(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_ref, 'factor', 1)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_ref(mod),3) ...
    ' Hz with global model'])

t2f = toc(t2);