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

%% SUBSTRUCTURING    

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

%At front of the heart
x_front = find(myMesh5.nodes(:,1)==1);
y_bc = find(abs(myMesh5.nodes(:,2))<=0.35);
z_bc = find(abs(myMesh5.nodes(:,3))<=0.35);

indices_front1 = intersect(x_front,y_bc);
indices_front = intersect(indices_front1,z_bc);

myMesh5.set_essential_boundary_condition(indices_front,1:3,0)

%At behind of the heart

x_behind = find(myMesh5.nodes(:,1)==0);

indices_behind1 = intersect(x_behind,y_bc);
indices_behind = intersect(indices_behind1,z_bc);

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
y_bc = find(abs(nodes_ref(:,2))<=0.35);
z_bc = find(abs(nodes_ref(:,3))<=0.35);

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

Assembly_ref.DATA.M = Assembly_ref.mass_matrix();
u0 = zeros( Mesh_ref.nDOFs, 1);
[Assembly_ref.DATA.K,~] = Assembly_ref.tangent_stiffness_and_force(u0);

Mc_ref = Assembly_ref.constrain_matrix(Assembly_ref.DATA.M);
Kc_ref = Assembly_ref.constrain_matrix(Assembly_ref.DATA.K);

figure
hold on
PlotMesh(nodes_ref, elements_ref, 0);
legend('Mesh ref')

%% VIBRATION MODES                                                

% EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

mod = 2;

%Substructuring

[~,om] = eigs(PrimalSub.DATA.Kc,PrimalSub.DATA.Mc, n_VMs, 'SM');
om = sort(sqrt(diag(om)));

%Global model

[~,om_ref] = eigs(Kc_ref, Mc_ref, n_VMs, 'SM');
om_ref = sort(sqrt(diag(om_ref)));

%% Reduction 

[M_hcb,K_hcb,V_hcb,L_hcb] = CraigBamptonReduction(PrimalSub,1000);

[~,om_hcb] = eigs(K_hcb,M_hcb, n_VMs, 'SM');
om_hcb = sort(diag(sqrt(om_hcb)));
%% Defining the force________________________________________________________

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

Fs = {f1,f2,f3,f4,f5};

F = L_to_global(PrimalSub,Fs);

%With global reference model

F_ref = zeros(Mesh_ref.nDOFs,1);

DOFs = Mesh_ref.get_DOF_from_location(loc1);
F_ref(DOFs(2)) = -S; %Force of Sub1

DOFs = Mesh_ref.get_DOF_from_location(loc2);
F_ref(DOFs(3)) = S; %Force of Sub2

DOFs = Mesh_ref.get_DOF_from_location(loc3);
F_ref(DOFs(2)) = S; %Force of Sub3

DOFs = Mesh_ref.get_DOF_from_location(loc4);
F_ref(DOFs(3)) = -S; %Force of Sub4

%% Time integration

%Constructing damping matrices

ksis = 0.1*ones(5,1);
A = [ones(5,1)./om/2 om.*ones(5,1)/2];
least_squares = (A'*A)\A'*ksis;
C = least_squares(1)*PrimalSub.DATA.M + least_squares(2)*PrimalSub.DATA.K;


A_ref = [ones(5,1)./om_ref/2 om_ref.*ones(5,1)/2];
least_squares_ref = (A_ref'*A_ref)\A_ref'*ksis;
Assembly_ref.DATA.C = least_squares_ref(1)*Assembly_ref.DATA.M + least_squares_ref(2)*Assembly_ref.DATA.K;


A_hcb = [ones(5,1)./om_hcb/2 om_hcb.*ones(5,1)/2];
least_squares_hcb = (A_hcb'*A_hcb)\A_hcb'*ksis;
C_hcb = least_squares_hcb(1)*M_hcb + least_squares_hcb(2)*K_hcb;
[nDOFhcb,~] = size(M_hcb);


omega_ext_ref = mean(om_ref(1:2));
omega_ext = mean(om(1:2));
omega_ext_hcb = mean(om_hcb(1:2));

T = 2*pi/omega_ext;
T_ref = 2*pi/omega_ext_ref;
T_hcb = 2*pi/omega_ext_hcb;

amplification_factor = 1;

%Forcing function
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);
F_ext_ref = @(t) amplification_factor * F_ref * sin(omega_ext_ref * t);


F_ext_hcb = applying_force_hcb(PrimalSub,Fs,V_hcb,L_hcb);
F_ext_hcb = -F_ext_hcb;
F_ext_hcb = @(t) amplification_factor * F_ext_hcb * sin(omega_ext_hcb * t);

% Initial condition: equilibrium
u0 = zeros(PrimalSub.nDOFglobal, 1);
v0 = zeros(PrimalSub.nDOFglobal, 1);
a0 = zeros(PrimalSub.nDOFglobal, 1); 

q0 = PrimalSub.constrain_vector(u0);
qd0 = PrimalSub.constrain_vector(v0);
qdd0 = PrimalSub.constrain_vector(a0); 

q0_ref = Assembly_ref.constrain_vector(u0);
qd0_ref = Assembly_ref.constrain_vector(v0);
qdd0_ref = Assembly_ref.constrain_vector(a0);

q0_hcb = zeros(nDOFhcb, 1);
qd0_hcb = zeros(nDOFhcb, 1);
qdd0_hcb = zeros(nDOFhcb, 1); 

% time step for integration
h = T/50;

TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
TI_lin_ref = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
TI_lin_hcb = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,PrimalSub,F_ext);
residual_lin_ref = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,Assembly_ref,F_ext_ref);
residual_HCB = @(q,qd,qdd,t)residual_HCB(q,qd,qdd, t, M_hcb,C_hcb, K_hcb, F_ext_hcb);

% Linearized Time Integration
tmax = 10*T; 
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);
TI_lin_ref.Integrate(q0_ref,qd0_ref,qdd0_ref,tmax,residual_lin_ref);
TI_lin_hcb.Integrate(q0_hcb,qd0_hcb,qdd0_hcb,tmax,residual_HCB);

% obtain full solution
TI_lin.Solution.u = PrimalSub.unconstrain_vector(TI_lin.Solution.q);
TI_lin_ref.Solution.u = Assembly_ref.unconstrain_vector(TI_lin_ref.Solution.q);

[~,s] = size(TI_lin_ref.Solution.u);

u_hcb_wrong = zeros(PrimalSub.nDOFglobal,s);
for i=1:s
    us_hcb_i = corr_reduction_indices(PrimalSub,V_hcb,L_hcb,TI_lin_hcb.Solution.q(:,i));
    
    u_hcb_wrong(:,i) = L_to_global(PrimalSub,us_hcb_i);
end

%% Differences




u = zeros(PrimalSub.nDOFglobal,s);
u_hcb = zeros(PrimalSub.nDOFglobal,s);
corr_gmsh_indices = corr_gmsh_indices(PrimalSub,Mesh_ref);

for i=2:s
    u(:,i) = reindex_vector(corr_gmsh_indices,TI_lin.Solution.u(:,i));
    u_hcb(:,i) = reindex_vector(corr_gmsh_indices,u_hcb_wrong(:,i));
end

u_norm = zeros(PrimalSub.nDOFglobal,s);
u_norm_ref = zeros(PrimalSub.nDOFglobal,s);
u_norm_hcb = zeros(PrimalSub.nDOFglobal,s);

angles = zeros(s,1);
angles_hcb = zeros(s,1);

for i=2:s
    u_norm(:,i) = u(:,i)/norm(u(:,i));
    u_norm_ref(:,i) = TI_lin_ref.Solution.u(:,i)/norm(TI_lin_ref.Solution.u(:,i));
    u_norm_hcb(:,i) = u_hcb(:,i)/norm(u_hcb(:,i));
    
    angles(i) = acos(u_norm(:,i)'*u_norm_ref(:,i))*180/pi;
    angles_hcb(i) = acos(u_norm_hcb(:,i)'*u_norm_ref(:,i))*180/pi;
end

mean_angle = mean(angles);
mean_angle_hcb = mean(angles_hcb);
