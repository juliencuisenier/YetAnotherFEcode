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

TotalMesh = load_gmsh2("RotorSketch.msh", -1); %reading the .msh file

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

myMesh5.set_essential_boundary_condition(Submeshes{5,3}{3},1:3,0)
myMesh5.set_essential_boundary_condition(Submeshes{5,3}{6},1:3,0)

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

Mesh_ref.set_essential_boundary_condition(nset_ref{3},1:3,0)
Mesh_ref.set_essential_boundary_condition(nset_ref{6},1:3,0)

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

[M_cb,K_cb,V_cb,L_cb] = CraigBamptonReduction(PrimalSub,1000);

[~,om_cb] = eigs(K_cb,M_cb, n_VMs, 'SM');
om_cb = sort(diag(sqrt(om_cb)));
%% Defining the force________________________________________________________

%With substructuring

S = 1e6;

%Sub1
loc1 = [-2.5 -0.15 0.5];
DOFs = PrimalSub.Substructures(1).Mesh.get_DOF_from_location(loc1);
f1=zeros(myMesh1.nDOFs,1);
f1(DOFs(2)) = S;

%Sub2
loc2 = [-0.15 2.5 0.5];
DOFs = PrimalSub.Substructures(2).Mesh.get_DOF_from_location(loc2);
f2=zeros(myMesh2.nDOFs,1);
f2(DOFs(1)) = S;

%Sub3
loc3 = [2.5 0.15 0.5];
DOFs = PrimalSub.Substructures(3).Mesh.get_DOF_from_location(loc3);
f3=zeros(myMesh3.nDOFs,1);
f3(DOFs(2)) = -S;

%Sub4
loc4 = [0.15 -2.5 0.5];
DOFs = PrimalSub.Substructures(4).Mesh.get_DOF_from_location(loc4);
f4=zeros(myMesh4.nDOFs,1);
f4(DOFs(1)) = -S;

%Sub5
f5=zeros(myMesh5.nDOFs,1);

Fs = {f1,f2,f3,f4,f5};

F = L_to_global(PrimalSub,Fs);

%With global reference model

F_ref = zeros(Mesh_ref.nDOFs,1);

DOFs = Mesh_ref.get_DOF_from_location(loc1);
F_ref(DOFs(2)) = S; %Force of Sub1

DOFs = Mesh_ref.get_DOF_from_location(loc2);
F_ref(DOFs(1)) = S; %Force of Sub2

DOFs = Mesh_ref.get_DOF_from_location(loc3);
F_ref(DOFs(2)) = -S; %Force of Sub3

DOFs = Mesh_ref.get_DOF_from_location(loc4);
F_ref(DOFs(1)) = -S; %Force of Sub4

%% Time integration

%Constructing damping matrices

ksis = 0.1*ones(5,1);
A = [ones(5,1)./om/2 om.*ones(5,1)/2];
least_squares = (A'*A)\A'*ksis;
PrimalSub.DATA.C = least_squares(1)*PrimalSub.DATA.M + least_squares(2)*PrimalSub.DATA.K;


A_ref = [ones(5,1)./om_ref/2 om_ref.*ones(5,1)/2];
least_squares_ref = (A_ref'*A_ref)\A_ref'*ksis;
Assembly_ref.DATA.C = least_squares_ref(1)*Assembly_ref.DATA.M + least_squares_ref(2)*Assembly_ref.DATA.K;


A_cb = [ones(5,1)./om_cb/2 om_cb.*ones(5,1)/2];
least_squares_cb = (A_cb'*A_cb)\A_cb'*ksis;
C_cb = least_squares_cb(1)*M_cb + least_squares_cb(2)*K_cb;
nDOFcb = size(M_cb,1);


omega_ext_ref = mean(om_ref(1:2));
omega_ext = mean(om(1:2));
omega_ext_cb = mean(om_cb(1:2));

T = 2*pi/omega_ext;
T_ref = 2*pi/omega_ext_ref;
T_cb = 2*pi/omega_ext_cb;

amplification_factor = 2;

%Forcing function
F_ext = @(t) amplification_factor * F * sin(omega_ext * t);
F_ext_ref = @(t) amplification_factor * F_ref * sin(omega_ext_ref * t);


F_ext_cb = applying_force_cb(PrimalSub,Fs,V_cb,L_cb);
F_ext_cb = @(t) amplification_factor * F_ext_cb * sin(omega_ext_cb * t);

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

q0_cb = zeros(nDOFcb, 1);
qd0_cb = zeros(nDOFcb, 1);
qdd0_cb = zeros(nDOFcb, 1); 

% time step for integration
h = T/5;

TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
TI_lin_ref = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
TI_lin_cb = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,PrimalSub,F_ext);
residual_lin_ref = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,Assembly_ref,F_ext_ref);
residual_cb = @(q,qd,qdd,t)residual_linear_cb(q,qd,qdd, t, M_cb,C_cb, K_cb, F_ext_cb);

% Linearized Time Integration
tmax = 10*T; 
TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);
TI_lin_ref.Integrate(q0_ref,qd0_ref,qdd0_ref,tmax,residual_lin_ref);
TI_lin_cb.Integrate(q0_cb,qd0_cb,qdd0_cb,tmax,residual_cb);

% obtain full solution
TI_lin.Solution.u = PrimalSub.unconstrain_vector(TI_lin.Solution.q);
TI_lin_ref.Solution.u = Assembly_ref.unconstrain_vector(TI_lin_ref.Solution.q);

[~,s] = size(TI_lin_ref.Solution.u);

u_cb_wrong = zeros(PrimalSub.nDOFglobal,s);
for i=1:s
    u_cb_wrong(:,i) = converter_reducted_vector(PrimalSub,V_cb,L_cb,TI_lin_cb.Solution.q(:,i));
end

%% Differences


u = zeros(PrimalSub.nDOFglobal,s);
u_cb = zeros(PrimalSub.nDOFglobal,s);
corr_gmsh_indices = corr_gmsh_indices(PrimalSub,Mesh_ref);

for i=2:s
    u(:,i) = reindex_vector(corr_gmsh_indices,TI_lin.Solution.u(:,i));
    u_cb(:,i) = reindex_vector(corr_gmsh_indices,u_cb_wrong(:,i));
end

u_norm = zeros(PrimalSub.nDOFglobal,s);
u_norm_ref = zeros(PrimalSub.nDOFglobal,s);
u_norm_cb = zeros(PrimalSub.nDOFglobal,s);

angles = zeros(s,1);
angles_cb = zeros(s,1);

for i=2:s
    u_norm(:,i) = u(:,i)/norm(u(:,i));
    u_norm_ref(:,i) = TI_lin_ref.Solution.u(:,i)/norm(TI_lin_ref.Solution.u(:,i));
    u_norm_cb(:,i) = u_cb(:,i)/norm(u_cb(:,i));
    
    angles(i) = acos(u_norm(:,i)'*u_norm_ref(:,i)*0.999999)*180/pi;
    angles_cb(i) = acos(u_norm_cb(:,i)'*u_norm_ref(:,i))*180/pi;
end

mean_angle = mean(angles);
mean_angle_cb = mean(angles_cb);

%% Plot
t = TI_lin.Solution.time;
dof = 223; %taking a random dof 
figure
plot(t,TI_lin_ref.Solution.u(dof,:),'k.-')
hold on
plot(t,u(dof,:),'b-')
plot(t,u_cb(dof,:),'r--')
legend('Reference model','Substructuring','Reduced model')
