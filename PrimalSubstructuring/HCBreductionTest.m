% Test of the PrimalSubstructuring class on a simili-rotor
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

x_front = find(Submeshes{5,1}(:,1)==1);
y_bc = find(abs(Submeshes{5,1}(:,2))<=0.5);
z_bc = find(abs(Submeshes{5,1}(:,3))<=0.5);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);


x_behind = find(Submeshes{5,1}(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

myMesh5.set_essential_boundary_condition(indices_front,1:3,0)
myMesh5.set_essential_boundary_condition(indices_behind,1:3,0)



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

%% REFERENCE MODEL_________________________________________________________
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

%% VIBRATION MODES                                                

% EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

mod = 2;

%Substructuring

[V0,om] = eigs(PrimalSub.DATA.Kc,PrimalSub.DATA.Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);

V0  = PrimalSub.unconstrain_vector(V0);
V0s = L_to_local(PrimalSub,V0);

figure
hold on
for jSub=1:PrimalSub.nSubs
    
    nodalDef = reshape(V0s{jSub}(:,mod),3,[]).';
    jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
    jElements = PrimalSub.Elements{jSub};
    PlotFieldonDeformedMesh(jMesh, jElements, nodalDef, 'factor', 10)
    
end
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ...
    ' Hz with substructuring'])


%Global model

[V0_ref,om_ref] = eigs(Kc_ref, Mc_ref, n_VMs, 'SM');
[f0_ref,ind_ref] = sort(sqrt(diag(om_ref))/2/pi);
V0_ref = V0_ref(:,ind_ref);

V0_ref = Assembly_ref.unconstrain_vector(V0_ref);

figure
hold on
nodalDef_ref = reshape(V0_ref(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_ref, 'factor', 10)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_ref(mod),3) ...
    ' Hz with global model'])

%% HCB reduction method

[M_hcb,K_hcb,L_hcb] = CraigBamptonReduction(PrimalSub,1000);

[V0_hcb,om] = eigs(M_hcb,M_hcb, n_VMs, 'SM');
[f0_hcb,ind] = sort(sqrt(diag(om))/2/pi);
V0_hcb = V0_hcb(:,ind);

V0_hcb  = PrimalSub.unconstrain_vector(V0_hcb);

Vs_hcb = {};

for iSub=1:PrimalSub.nSubs
    Vs_hcb{iSub} = L_hcb{iSub}*V0_hcb;

    figure
    hold on

    
    nodalDef = reshape(V0s{iSub}(:,mod),3,[]).';
    iMesh = PrimalSub.Substructures(iSub).Mesh.nodes;
    iElements = PrimalSub.Elements{iSub};
    PlotFieldonDeformedMesh(iMesh, iElements, nodalDef, 'factor', 10)
    
end
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_hcb(mod),3) ...
    ' Hz with HCB reduction'])

