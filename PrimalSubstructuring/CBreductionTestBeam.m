% Test of the PrimalSubstructuring class on a clambed-clambed beam
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

TotalMesh = load_gmsh2("geomBeam.msh", -1); %reading the .msh file

Submeshes = Submeshing(TotalMesh);


%Mesh 1___________________________________________________________________

myMesh1 = Mesh(Submeshes{1,1});
myMesh1.create_elements_table(Submeshes{1,2}, myElementConstructor);
myMesh1.set_essential_boundary_condition(Submeshes{1,3}{5},1:3,0)

%Mesh 2____________________________________________________________________

myMesh2 = Mesh(Submeshes{2});
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

globalIndices = get_globalIndices(Submeshes);
 
t1 = tic;
PrimalSub = PrimalSubstructuring([Assembly1 Assembly2],globalIndices,[]);
t_sub_Assembly = toc(t1);

%% REFERENCE MODEL

[nodes_ref,elements_ref,nset_ref] = extract_gmsh(TotalMesh);

Mesh_ref = Mesh(nodes_ref);
Mesh_ref.create_elements_table(elements_ref, myElementConstructor)
Mesh_ref.set_essential_boundary_condition(nset_ref{2},1:3,0)
Mesh_ref.set_essential_boundary_condition(nset_ref{5},1:3,0)

t1=tic;
Assembly_ref = Assembly(Mesh_ref);

M_ref = Assembly_ref.mass_matrix();
u0 = zeros( Mesh_ref.nDOFs, 1);
[K_ref,~] = Assembly_ref.tangent_stiffness_and_force(u0);

Mc_ref = Assembly_ref.constrain_matrix(M_ref);
Kc_ref = Assembly_ref.constrain_matrix(K_ref);

Assembly_ref.DATA.M = M_ref;
Assembly_ref.DATA.Mc = Mc_ref;
Assembly_ref.DATA.K = K_ref;
Assembly_ref.DATA.Kc = Kc_ref;

t_ref_Assembly = toc(t1);

%% CB Reduction method
t1 = tic;
[M_cb,K_cb,V_cb,L_cb] = CraigBamptonReduction(PrimalSub,3000);
t_cb_Assembly = toc(t1);

%% VIBRATION MODES                                                  

% % EIGENMODES ______________________________________________________________
n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
% 
% 
% t_sub = tic;
% 
% [V0_sub,om] = eigs(PrimalSub.DATA.Kc,PrimalSub.DATA.Mc, n_VMs, 'SM');
% [f0,ind] = sort(sqrt(diag(om))/2/pi);
% V0_sub = V0_sub(:,ind);
% 
% t_sub_VM = toc(t_sub);
% 
% V0_sub  = PrimalSub.unconstrain_vector(V0_sub);
% V0s = L_to_local(PrimalSub,V0_sub);
% 
% figure
% hold on
% for jSub=1:PrimalSub.nSubs
%     
%     nodalDef = reshape(V0s{jSub}(:,mod),3,[]).';
%     jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
%     jElements = PrimalSub.Elements{jSub};
%     PlotFieldonDeformedMesh(jMesh, jElements, nodalDef, 'factor', 1)
%     
% end
% colormap jet
% title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ...
%     ' Hz with substructuring'])


%Global model

t1 = tic;

[V0_ref,om_ref] = eigs(Kc_ref, Mc_ref, n_VMs, 'SM');
[f0_ref,ind_ref] = sort(sqrt(diag(om_ref))/2/pi);
V0_ref = V0_ref(:,ind_ref);

t_ref_VM = toc(t1);

V0_ref = Assembly_ref.unconstrain_vector(V0_ref);




% Craig-Bampton reduction

t1 = tic;

[V0_cb,om] = eigs(K_cb,M_cb, n_VMs, 'SM');
[f0_cb,ind] = sort(sqrt(diag(om))/2/pi);
V0_cb = V0_cb(:,ind);

V0_cb_converted = zeros(PrimalSub.nDOFglobal,n_VMs);

for i=1:n_VMs
    V0_cb_converted(:,i) = converter_reducted_vector(PrimalSub,V_cb,L_cb,V0_cb(:,i));
end

t_cb_VM =  toc(t1);

% %corr_gmsh_indices = corr_gmsh_indices(PrimalSub,Mesh_ref);
% V0_cb_converted(:,2) = reindex_vector(corr_gmsh_indices,V0_cb_converted(:,2));
% 
% figure
% hold on
% nodalDef_hcb = reshape(V0_cb_converted(:,2),3,[]).';
% PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_hcb, 'factor', 1)
% colormap jet
% title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_cb(mod),3) ...
%     ' Hz with HCB reduction'])

%% Defining the force
S = 1e6;

loc = [0.1 0.5 0.3];

% Substructuring
f1 = zeros(myMesh1.nDOFs,1);
DOFs = myMesh1.get_DOF_from_location(loc);
f1(DOFs(3)) = S;

f2 = zeros(myMesh2.nDOFs,1);

Fext_sub = {f1,f2};

F_sub = L_to_global(PrimalSub,Fext_sub);

% Global model

Fext_ref = zeros(Mesh_ref.nDOFs,1);
DOFs = Mesh_ref.get_DOF_from_location(loc);
Fext_ref(DOFs(3)) = S;

%% Static resolution

% t1 = tic;
% u_sub = static_resolution(PrimalSub,F_sub);
% t_sub_static = toc(t1);

% Global model

t1 = tic;
Fextc_ref = Assembly_ref.constrain_vector(Fext_ref);

uc_ref = Kc_ref\Fextc_ref;
u_ref = Assembly_ref.unconstrain_vector(uc_ref);

t_ref_static = toc(t1);

% CB reduction



t1 = tic;
Fext_cb = applying_force_cb(PrimalSub,Fext_sub,V_cb,L_cb);
u_cb = K_cb\Fext_cb;
u_cb = converter_reducted_vector(PrimalSub,V_cb,L_cb,u_cb);

t_cb_static = toc(t1);

%% Damping matrices

om = 2*pi*f0_ref;
om_ref = 2*pi*f0_ref;
om_cb = 2*pi*f0_cb;

ksis = 0.1*ones(n_VMs,1);
A = [ones(5,1)./om/2 om.*ones(n_VMs,1)/2];
least_squares = (A'*A)\A'*ksis;
PrimalSub.DATA.C = least_squares(1)*PrimalSub.DATA.M + least_squares(2)*PrimalSub.DATA.K;
PrimalSub.DATA.Cc = PrimalSub.constrain_matrix(PrimalSub.DATA.C);


A_ref = [ones(n_VMs,1)./om_ref/2 om_ref.*ones(n_VMs,1)/2];
least_squares_ref = (A_ref'*A_ref)\A_ref'*ksis;
Assembly_ref.DATA.C = least_squares_ref(1)*Assembly_ref.DATA.M + least_squares_ref(2)*Assembly_ref.DATA.K;
Assembly_ref.DATA.Cc = Assembly_ref.constrain_matrix(Assembly_ref.DATA.C);

A_cb = [ones(n_VMs,1)./om_cb/2 om_cb.*ones(n_VMs,1)/2];
least_squares_cb = (A_cb'*A_cb)\A_cb'*ksis;
C_cb = least_squares_cb(1)*M_cb + least_squares_cb(2)*K_cb;
nDOFcb = size(M_cb,1);

%% FRFs

fmin = 0;
fmax = 2500;
fn = 100;

% Substructuring
% t1 = tic;
% FRF_sub = frf_substructuring(PrimalSub,Fext_sub,DOFs,fmin,fmax,fn);
% t_FRF_sub = toc(t1);

% Full model
t1 = tic;

FRF_ref = zeros(Mesh_ref.nDOFs,fn);
freq = fmin:(fmax-fmin)/(fn-1):fmax;
om = 2*pi*freq;

for i=1:fn
    G = -om(i)^2*Assembly_ref.DATA.Mc +1i*om(i)*Assembly_ref.DATA.Cc +Assembly_ref.DATA.Kc;
    FRF_ref(:,i) = Assembly_ref.unconstrain_vector(G\Fextc_ref);
end
t_FRF_ref = toc(t1);

t1 = tic;
FRF_cb = frf_cb(PrimalSub,M_cb,C_cb,K_cb,V_cb,L_cb,Fext_sub,DOFs,fmin,fmax,fn);
t_FRF_cb = toc(t1);

%% Time integration

omega_ext_ref = mean(f0_ref(1:2)*2*pi);
omega_ext = mean(f0_ref(1:2)*2*pi);
omega_ext_cb = mean(f0_ref(1:2)*2*pi);

T = 2*pi/omega_ext;

amplification_factor = 2;

%Forcing function
F_ext_sub = @(t) amplification_factor * F_sub * sin(omega_ext * t);
F_ext_ref = @(t) amplification_factor * Fext_ref * sin(omega_ext_ref * t);


F_ext_cb = applying_force_cb(PrimalSub,Fext_sub,V_cb,L_cb);
F_ext_cb = @(t) amplification_factor * F_ext_cb * sin(omega_ext_ref * t);

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

% TI_lin = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
TI_lin_ref = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);
TI_lin_cb = ImplicitNewmark('timestep',h,'alpha',0.005,'linear',true);

% Linear Residual evaluation function handle
% residual_lin = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,PrimalSub,F_ext_sub);
residual_lin_ref = @(q,qd,qdd,t)residual_linear(q,qd,qdd,t,Assembly_ref,F_ext_ref);
residual_lin_cb = @(q,qd,qdd,t)residual_linear_cb(q,qd,qdd, t, M_cb,C_cb, K_cb, F_ext_cb);

% Linearized Time Integration
tmax = 10*T;

% t1 = tic;
% TI_lin.Integrate(q0,qd0,qdd0,tmax,residual_lin);
% t_TI_sub = toc(t1);

t1 = tic;
TI_lin_ref.Integrate(q0_ref,qd0_ref,qdd0_ref,tmax,residual_lin_ref);
t_TI_ref = toc(t1);

t1 = tic;
TI_lin_cb.Integrate(q0_cb,qd0_cb,qdd0_cb,tmax,residual_lin_cb);
t_TI_cb = toc(t1);

% obtain full solution
% TI_lin.Solution.u = PrimalSub.unconstrain_vector(TI_lin.Solution.q);
TI_lin_ref.Solution.u = Assembly_ref.unconstrain_vector(TI_lin_ref.Solution.q);

s = size(TI_lin_ref.Solution.u,2);

u_TI_cb = zeros(PrimalSub.nDOFglobal,s);
for i=1:s
    u_TI_cb(:,i) = converter_reducted_vector(PrimalSub,V_cb,L_cb,TI_lin_cb.Solution.q(:,i));
end

%% %%%%%%%%%%%%%%%%%%%% Differences %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr_gmsh_indices = corr_gmsh_indices(PrimalSub,Mesh_ref);

%% VMs

% angle_VM_ref_sub = zeros(1,n_VMs);
angle_VM_ref_cb = zeros(1,n_VMs);

for i=1:n_VMs
%     V0_sub(:,i) = reindex_vector(corr_gmsh_indices,V0_sub(:,i));
    V0_cb_converted(:,i) = reindex_vector(corr_gmsh_indices,V0_cb_converted(:,i));
    
%     normV0_sub = V0_sub(:,i)/norm(V0_sub(:,i));
    normV0_ref = V0_ref(:,i)/norm(V0_ref(:,i));
    normV0_cb = V0_cb_converted(:,i)/norm(V0_cb_converted(:,i));
    
%     angle_VM_ref_sub(:,i) = acos(normV0_sub'*normV0_ref)*180/pi;
    angle_VM_ref_cb(:,i) = acos(normV0_cb'*normV0_ref)*180/pi;  
end


% mean_angle_VM_ref_sub = mean(angle_VM_ref_sub);
mean_angle_VM_ref_cb = mean(angle_VM_ref_cb);

%% Static resolution

% u_sub = reindex_vector(corr_gmsh_indices,u_sub);
u_cb = reindex_vector(corr_gmsh_indices,u_cb);

% norm_u_sub = u_sub/norm(u_sub);
norm_u_ref = u_ref/norm(u_ref);
norm_u_cb = u_cb/norm(u_cb);

% angle_u_ref_sub = acos(norm_u_sub'*norm_u_ref)*180/pi;
angle_u_ref_cb = acos(norm_u_cb'*norm_u_ref)*180/pi;  

%% FRFs

% angle_FRF_ref_sub = zeros(1,fn);
angle_FRF_ref_cb = zeros(1,fn);



for i=1:fn
%     FRF_sub(:,i) = reindex_vector(corr_gmsh_indices,abs(FRF_sub(:,i)));
    FRF_cb(:,i) = reindex_vector(corr_gmsh_indices,abs(FRF_cb(:,i)));
    
%     normFRF_sub = FRF_sub(:,i)/norm(FRF_sub(:,i));
    normFRF_ref = abs(FRF_ref(:,i))/norm(abs(FRF_ref(:,i)));
    normFRF_cb = FRF_cb(:,i)/norm(FRF_cb(:,i));
    
   
%     angle_FRF_ref_sub(:,i) = acos(normFRF_sub'*normFRF_ref*0.99999999)*180/pi;
    angle_FRF_ref_cb(:,i) = acos(normFRF_cb'*normFRF_ref)*180/pi;  
end


% mean_angle_FRF_ref_sub = mean(angle_FRF_ref_sub);
mean_angle_FRF_ref_cb = mean(angle_FRF_ref_cb);

%% Time Integration

% Reminder : s is the number of time step, i.e the number of columns of the
% TI solution
% angle_TI_ref_sub = zeros(1,s);
angle_TI_ref_cb = zeros(1,s);

% u = zeros(PrimalSub.nDOFglobal,s);

for i=2:s
%     u(:,i) = reindex_vector(corr_gmsh_indices,TI_lin.Solution.u(:,i));
    u_TI_cb(:,i) = reindex_vector(corr_gmsh_indices,u_TI_cb(:,i));
    
%     norm_u_TI_sub = u(:,i)/norm(u(:,i));
    norm_u_TI_ref = TI_lin_ref.Solution.u(:,i)/norm(TI_lin_ref.Solution.u(:,i));
    norm_u_TI_cb = u_TI_cb(:,i)/norm(u_TI_cb(:,i));
    
%     angle_TI_ref_sub(i) = acos(norm_u_TI_sub'*norm_u_TI_ref*0.999999)*180/pi;
    angle_TI_ref_cb(i) = acos(norm_u_TI_cb'*norm_u_TI_ref)*180/pi;
end

% mean_angle_TI_ref_sub = mean(angle_TI_ref_sub);
mean_angle_TI_ref_cb = mean(angle_TI_ref_cb);
