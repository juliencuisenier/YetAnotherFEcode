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


TotalMesh = load_gmsh2("BeamMass2.msh", -1); %reading the .msh file


[nodes_ref,elements_ref,nset_ref] = extract_gmsh(TotalMesh);

Mesh_ref = Mesh(nodes_ref);
Mesh_ref.create_elements_table(elements_ref, myElementConstructor)

%Boundaries conditions

Mesh_ref.set_essential_boundary_condition(nset_ref{1},1:3,0)


%Assembly

Assembly_ref = Assembly(Mesh_ref);

M_ref = Assembly_ref.mass_matrix();
u0 = zeros( Mesh_ref.nDOFs, 1);
[K_ref,~] = Assembly_ref.tangent_stiffness_and_force(u0);

Mc_ref = Assembly_ref.constrain_matrix(M_ref);
Kc_ref = Assembly_ref.constrain_matrix(K_ref);

Assembly_ref.DATA.M = M_ref;
Assembly_ref.DATA.Mc = Assembly_ref.constrain_matrix(M_ref);
Assembly_ref.DATA.K = K_ref;
Assembly_ref.DATA.Kc = Assembly_ref.constrain_matrix(K_ref);

figure
hold on
PlotMesh(nodes_ref, elements_ref, 0);
legend('Mesh ref')


%% VIBRATION MODES                                                

% EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 

mod = 2;

[V0_ref,om_ref] = eigs(Kc_ref, Mc_ref, n_VMs, 'SM');
[f0_ref,ind_ref] = sort(sqrt(diag(om_ref))/2/pi);
om_ref = sort(sqrt(diag(om_ref)));
V0_ref = V0_ref(:,ind_ref);

V0_ref = Assembly_ref.unconstrain_vector(V0_ref);

figure
hold on
nodalDef_ref = reshape(V0_ref(:,mod),3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_ref, 'factor', 10)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_ref(mod),3) ...
    ' Hz with global model'])

%% Constructing damping matrices

ksis = 0.1*ones(5,1);


A_ref = [ones(5,1)./om_ref/2 om_ref.*ones(5,1)/2];
least_squares_ref = (A_ref'*A_ref)\A_ref'*ksis;
Assembly_ref.DATA.C = least_squares_ref(1)*Assembly_ref.DATA.M + least_squares_ref(2)*Assembly_ref.DATA.K;
Assembly_ref.DATA.Cc = Assembly_ref.constrain_matrix(Assembly_ref.DATA.C);

%% Defining the excitation force___________________________________________

S = 1e6;

loc = [5.5 0 0.5];
DOFs = Assembly_ref.Mesh.get_DOF_from_location(loc);
F_ref=zeros(Mesh_ref.nDOFs,1);
F_ref(DOFs(3)) = -S;

Fc_ref = Assembly_ref.constrain_vector(F_ref);


%% FRF

DOFs_ref = Assembly_ref.Mesh.get_DOF_from_location([5.5 0.5 0]);

fmax = 100;
fn = 100;

%Calculating the FRF og the global model
FRF_ref = zeros(Mesh_ref.nDOFs,fn);
freqs = 0:fmax/(fn-1):fmax;
omegas = 2*pi*freqs;
for i=1:fn
    G = -omegas(i)^2*Assembly_ref.DATA.Mc +1i*omegas(i)*Assembly_ref.DATA.Cc +Assembly_ref.DATA.Kc;
    FRF_ref(:,i) = Assembly_ref.unconstrain_vector(G\Fc_ref);
end


for xi=DOFs_ref(2)
   figure
   hold on
   subplot(2,1,1);
   plot(freqs,abs(FRF_ref(xi,:)));
   title(strcat('Module of the DOF ',num2str(xi),"'s FRF (global model)"));
   xlabel('Frequency');
   set(gca,'yscale','log');
   subplot(2,1,2);
   plot(freqs,unwrap(angle(FRF_ref(xi,:))))
   title(strcat('Phase of the DOF ',num2str(xi),"'s FRF [rad]"));
   xlabel('Frequency');
end

%% Time integration

omega_ext_ref = 10*2*pi;

T_ref = 2*pi/omega_ext_ref;

amplification_factor = 2;

%Forcing function
F_ext_ref = @(t) amplification_factor * F_ref * sin(omega_ext_ref * t);



% Initial condition: equilibrium
u0 = zeros(Mesh_ref.nDOFs, 1);
v0 = zeros(Mesh_ref.nDOFs, 1);
a0 = zeros(Mesh_ref.nDOFs, 1);

q0_ref = Assembly_ref.constrain_vector(u0);
qd0_ref = Assembly_ref.constrain_vector(v0);
qdd0_ref = Assembly_ref.constrain_vector(a0);


% time step for integration
h = T_ref/5;

TI_lin_ref = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin_ref = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,Assembly_ref,F_ext_ref);

% Linearized Time Integration
tmax = 10*T_ref; 
TI_lin_ref.Integrate(q0_ref,qd0_ref,qdd0_ref,tmax,residual_lin_ref);


% obtain full solution
TI_lin_ref.Solution.u = Assembly_ref.unconstrain_vector(TI_lin_ref.Solution.q);

t = TI_lin_ref.Solution.time;
figure
plot(t,TI_lin_ref.Solution.u(DOFs_ref(2),:),'k.-')
hold on


%% Rotor forces
Omega = [147 0 0]';

%For ref model

G_ref = Assembly_ref.coriolis_matrix(Omega);
Ksp_ref = Assembly_ref.spin_softening_matrix(Omega);

Assembly_ref.DATA.C = Assembly_ref.DATA.C + G_ref;
Assembly_ref.DATA.K = Assembly_ref.DATA.K - Ksp_ref;
Assembly_ref.DATA.Cc = Assembly_ref.constrain_matrix(Assembly_ref.DATA.C);
Assembly_ref.DATA.Kc = Assembly_ref.constrain_matrix(Assembly_ref.DATA.K);



%% VM with rotor forces

[V0_ref_r,om_ref_r] = eigs(Assembly_ref.DATA.Kc,Assembly_ref.DATA.Cc, n_VMs, 'SM');
[f0_ref_r,ind_ref] = sort(sqrt(diag(om_ref_r))/2/pi);
om_ref_r = sort(sqrt(diag(om_ref_r)));
V0_ref_r = V0_ref_r(:,ind_ref);

V0_ref_r = Assembly_ref.unconstrain_vector(V0_ref_r);

figure
hold on
nodalDef_ref = reshape(abs(V0_ref_r(:,mod)),3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_ref, 'factor', 10)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_ref_r(mod),3) ...
    ' Hz with global model'])


%% FRF rotor

DOFs_ref = Assembly_ref.Mesh.get_DOF_from_location([5.5 0.5 0]);

fmax = 100;
fn = 100;

%Calculating the FRF og the global model
FRF_ref = zeros(Mesh_ref.nDOFs,fn);
freqs = 0:fmax/(fn-1):fmax;
omegas = 2*pi*freqs;
for i=1:fn
    G = -omegas(i)^2*Assembly_ref.DATA.Mc +1i*omegas(i)*Assembly_ref.DATA.Cc +Assembly_ref.DATA.Kc;
    FRF_ref(:,i) = Assembly_ref.unconstrain_vector(G\Fc_ref);
end


for xi=DOFs_ref(2)
   figure
   hold on
   subplot(2,1,1);
   plot(freqs,abs(FRF_ref(xi,:)));
   title(strcat('Module of the DOF ',num2str(xi),"'s FRF (global model)"));
   xlabel('Frequency');
   set(gca,'yscale','log');
   subplot(2,1,2);
   plot(freqs,unwrap(angle(FRF_ref(xi,:))))
   title(strcat('Phase of the DOF ',num2str(xi),"'s FRF [rad]"));
   xlabel('Frequency');
end

%% Time integration

omega_ext_ref = 10*2*pi;

T_ref = 2*pi/omega_ext_ref;

amplification_factor = 2;

%Forcing function
F_ext_ref = @(t) amplification_factor * F_ref * sin(omega_ext_ref * t);



% Initial condition: equilibrium
u0 = zeros(Mesh_ref.nDOFs, 1);
v0 = zeros(Mesh_ref.nDOFs, 1);
a0 = zeros(Mesh_ref.nDOFs, 1);

q0_ref = Assembly_ref.constrain_vector(u0);
qd0_ref = Assembly_ref.constrain_vector(v0);
qdd0_ref = Assembly_ref.constrain_vector(a0);


% time step for integration
h = T_ref/5;

TI_lin_ref = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
residual_lin_ref = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,Assembly_ref,F_ext_ref);

% Linearized Time Integration
tmax = 10*T_ref; 
TI_lin_ref.Integrate(q0_ref,qd0_ref,qdd0_ref,tmax,residual_lin_ref);


% obtain full solution
TI_lin_ref.Solution.u = Assembly_ref.unconstrain_vector(TI_lin_ref.Solution.q);

t = TI_lin_ref.Solution.time;
figure
plot(t,TI_lin_ref.Solution.u(DOFs_ref(2),:),'k.-')
hold on


