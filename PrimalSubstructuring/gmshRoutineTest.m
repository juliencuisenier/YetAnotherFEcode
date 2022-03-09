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
myMesh1.set_essential_boundary_condition(Submeshes{1,4}{1},1:3,0)

figure
hold on
PlotMesh(Submeshes{1,1}, Submeshes{1,2}, 0);
legend('Mesh SS#1')

%Mesh 2____________________________________________________________________

myMesh2 = Mesh(Submeshes{2,1});
myMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);
myMesh2.set_essential_boundary_condition(Submeshes{2,4}{1},1:3,0)

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


 
% Test of the static resolution ___________________________________________

% Some arbitrary external forces, fextN corresponding to the external force
% applied on SubN


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