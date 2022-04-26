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

x_front = find(Submeshes{1,1}(:,1)==1);
y_bc = find(abs(Submeshes{1,1}(:,2))<=0.5);
z_bc = find(abs(Submeshes{1,1}(:,3))<=0.5);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);

myMesh1.set_essential_boundary_condition(indices_front,1:3,0)

x_behind = find(Submeshes{1,1}(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

myMesh1.set_essential_boundary_condition(indices_behind,1:3,0)

%Mesh 2 (Top blade)________________________________________________________

myMesh2 = Mesh(Submeshes{2,1});
myMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);

x_front = find(Submeshes{2,1}(:,1)==1);
y_bc = find(abs(Submeshes{2,1}(:,2))<=0.5);
z_bc = find(abs(Submeshes{2,1}(:,3))<=0.5);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);

myMesh2.set_essential_boundary_condition(indices_front,1:3,0)

x_behind = find(Submeshes{2,1}(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

myMesh2.set_essential_boundary_condition(indices_behind,1:3,0)

%Mesh 3 (Right blade)______________________________________________________

myMesh3 = Mesh(Submeshes{3,1});
myMesh3.create_elements_table(Submeshes{3,2}, myElementConstructor);

x_front = find(Submeshes{3,1}(:,1)==1);
y_bc = find(abs(Submeshes{3,1}(:,2))<=0.5);
z_bc = find(abs(Submeshes{3,1}(:,3))<=0.5);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);

myMesh3.set_essential_boundary_condition(indices_front,1:3,0)

x_behind = find(Submeshes{3,1}(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

myMesh3.set_essential_boundary_condition(indices_behind,1:3,0)

%Mesh 4 (Bottom blade)_____________________________________________________

myMesh4 = Mesh(Submeshes{4,1});
myMesh4.create_elements_table(Submeshes{4,2}, myElementConstructor);

x_front = find(Submeshes{4,1}(:,1)==1);
y_bc = find(abs(Submeshes{4,1}(:,2))<=0.5);
z_bc = find(abs(Submeshes{4,1}(:,3))<=0.5);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);

myMesh4.set_essential_boundary_condition(indices_front,1:3,0)

x_behind = find(Submeshes{4,1}(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

myMesh4.set_essential_boundary_condition(indices_behind,1:3,0)

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

%% STATIC RESOLUTION

%Defining the force________________________________________________________

%With substructuring

S = 1e6;

%Sub1
loc1 = [0.5 0.15 2.5];
DOFs = PrimalSub.Substructures(1).Mesh.get_DOF_from_location(loc1);
f1=zeros(myMesh1.nDOFs,1);
f1(DOFs(2)) = -S; %Force applied on the 2nd DOF of the 33rd node

%Sub2
loc2 = [0.5 2.5 -0.15];
DOFs = PrimalSub.Substructures(2).Mesh.get_DOF_from_location(loc2);
f2=zeros(myMesh2.nDOFs,1);
f2(DOFs(3)) = S; %Force applied on the 3th DOF of the 32nd node

%Sub3
loc3 = [0.5 -0.15 -2.5];
DOFs = PrimalSub.Substructures(3).Mesh.get_DOF_from_location(loc3);
f3=zeros(myMesh3.nDOFs,1);
f3(DOFs(2)) = S; %Force applied on the 2nd DOF of the 22nd node

%Sub4
loc4 = [0.5 -2.5 0.15];
DOFs = PrimalSub.Substructures(4).Mesh.get_DOF_from_location(loc4);
f4=zeros(myMesh4.nDOFs,1);
f4(DOFs(3)) = -S; %Force applied on the 3rd DOF of the 12th node

%Sub5
f5=zeros(myMesh5.nDOFs,1);

Fext = {f1,f2,f3,f4,f5};

%With global reference model

Fext_ref = zeros(Mesh_ref.nDOFs,1);

DOFs = Mesh_ref.get_DOF_from_location(loc1);
Fext_ref(DOFs(2)) = -S; %Force of Sub1

DOFs = Mesh_ref.get_DOF_from_location(loc2);
Fext_ref(DOFs(3)) = S; %Force of Sub2

DOFs = Mesh_ref.get_DOF_from_location(loc3);
Fext_ref(DOFs(2)) = S; %Force of Sub3

DOFs = Mesh_ref.get_DOF_from_location(loc4);
Fext_ref(DOFs(3)) = -S; %Force of Sub4

%Static resolutions________________________________________________________

%Using substructuring
u_sub = static_resolution(PrimalSub,Fext);

us = L_to_local(PrimalSub,u_sub);

figure
hold on
for jSub=1:PrimalSub.nSubs
    
    staticDef = reshape(us{jSub},3,[]).';
    jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
    jElements = PrimalSub.Elements{jSub};
    PlotFieldonDeformedMesh(jMesh, jElements, staticDef, 'factor', 10)
    
end
colormap jet
title('Deformation of the structure using substructuring')

%Using global model

Fextc_ref = Assembly_ref.constrain_vector(Fext_ref);

uc_ref = Kc_ref\Fextc_ref;

u_ref = Assembly_ref.unconstrain_vector(uc_ref);

figure
hold on
staticDef_ref = reshape(u_ref,3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, staticDef_ref, 'factor', 10)
colormap jet
title('Deformation of the structure using global model')

%% DEFECT MODES

%Using substructuring


[~,new_K] = PrimalSub.global_mass_stiffness(us);
new_Kc = PrimalSub.constrain_matrix(new_K);

[new_V0,new_om] = eigs(new_Kc, PrimalSub.DATA.Mc, n_VMs, 'SM');
[new_f0,new_ind] = sort(sqrt(diag(new_om))/2/pi);
new_V0 = new_V0(:,new_ind);

new_V0  = PrimalSub.unconstrain_vector(new_V0);
new_V0s = L_to_local(PrimalSub,new_V0);

figure
hold on
for jSub=1:PrimalSub.nSubs
    
    new_nodalDef = reshape(new_V0s{jSub}(:,mod),3,[]).';
    jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
    jElements = PrimalSub.Elements{jSub};
    PlotFieldonDeformedMesh(jMesh, jElements, new_nodalDef, 'factor', 10)
    
end
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(new_f0(mod),3) ...
    ' Hz using substructuring with defects'])

%With global model

[new_K_ref,~] = Assembly_ref.tangent_stiffness_and_force(u_ref);

new_Kc_ref = Assembly_ref.constrain_matrix(new_K_ref);

[new_V0_ref,new_om_ref] = eigs(new_Kc_ref, Mc_ref, n_VMs, 'SM');
[new_f0_ref,new_ind_ref] = sort(sqrt(diag(new_om_ref))/2/pi);
new_V0_ref = new_V0_ref(:,new_ind_ref);

new_V0_ref = Assembly_ref.unconstrain_vector(new_V0_ref);

figure
hold on
new_nodalDef_ref = reshape(new_V0_ref(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, new_nodalDef_ref, 'factor', 10)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(new_f0_ref(mod),3) ...
    ' Hz using global model with defects'])


%% Differences

corr_gmsh_indices = corr_gmsh_indices(PrimalSub,Mesh_ref);
V0 = reindex_vector(corr_gmsh_indices,V0(:,2));
u_sub = reindex_vector(corr_gmsh_indices,u_sub);
new_V0 = reindex_vector(corr_gmsh_indices,new_V0(:,2));

%Vibration modes
V0_norm = V0/norm(V0);
V0_ref_norm = V0_ref(:,2)/norm(V0_ref(:,2));

VM_deflects_angle = acos(V0_norm'*V0_ref_norm)*180/pi;

%Static deflections
u_sub_norm = u_sub/norm(u_sub);
u_ref_norm = u_ref/norm(u_ref);

static_deflects_angle = acos(u_sub_norm'*u_ref_norm)*180/pi;

%VM with defects
new_V0_norm = new_V0/norm(new_V0);
new_V0_ref_norm = new_V0_ref(:,2)/norm(new_V0_ref(:,2));

new_VM_deflects_angle = acos(new_V0_norm'*new_V0_ref_norm)*180/pi;
