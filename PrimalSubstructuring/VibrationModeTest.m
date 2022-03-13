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

TotalMesh = load_gmsh2("GeomBeam.msh", -1); %reading the .msh file

Submeshes = Submeshing(TotalMesh,2);


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


% EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

[V0,om] = eigs(PrimalSub.DATA.Kc, PrimalSub.DATA.Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);

 
V0  = PrimalSub.unconstrain_vector(V0);
%V0  = ROBnormalization( V0, 'displacement' );
V0s = deformation2substructs(PrimalSub,V0);   % Localize V0 to substructure s

% PLOT of the two substructures ___________________________________________
mod = 2;
figure
hold on
%PlotMesh(Submeshes{1,1}, Submeshes{1,2}, 0);
%PlotMesh(Submeshes{2,1}, Submeshes{2,2}, 0);
nodalDef1 = reshape(V0s{1}(:,mod),3,[]).';
nodalDef2 = reshape(V0s{2}(:,mod),3,[]).';
PlotFieldonDeformedMesh(Submeshes{1,1}, Submeshes{1,2}, nodalDef1, 'factor', 1)
PlotFieldonDeformedMesh(Submeshes{2,1}, Submeshes{2,2}, nodalDef2, 'factor', 1)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ...
    ' Hz with substructuring'])

% REFERENCE MODEL 

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
[f0_ref,ind_ref] = sort(sqrt(diag(om))/2/pi);
V0_ref = V0_ref(:,ind);
f0_full = f0;

V0_ref = Assembly_ref.unconstrain_vector(V0_ref);

figure
hold on
nodalDef_ref = reshape(V0_ref(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_ref, 'factor', 1)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_ref(mod),3) ...
    ' Hz with global model'])