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

z_front = find(Submeshes{5,1}(:,3)==1);
x_bc = find(abs(Submeshes{5,1}(:,1))<=0.5);
y_bc = find(abs(Submeshes{5,1}(:,2))<=0.5);

indices_front1 = intersect(z_front,x_bc);
indices_front = intersect(indices_front1,y_bc);


z_behind = find(Submeshes{5,1}(:,3)==0);

indices_behind1 = intersect(z_behind,x_bc);
indices_behind = intersect(indices_behind1,y_bc);

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

%% REFERENCE MODEL
[nodes_ref,elements_ref,nset_ref] = extract_gmsh(TotalMesh);

Mesh_ref = Mesh(nodes_ref);
Mesh_ref.create_elements_table(elements_ref, myElementConstructor)

%Boundaries conditions

%At front of the heart
z_front = find(nodes_ref(:,3)==1);
x_bc = find(abs(nodes_ref(:,1))<=0.5);
y_bc = find(abs(nodes_ref(:,2))<=0.5);

indices_front1 = intersect(z_front,x_bc);
indices_front = intersect(indices_front1,y_bc);

Mesh_ref.set_essential_boundary_condition(indices_front,1:3,0)

%At behind of the heart

z_behind = find(nodes_ref(:,3)==0);

indices_behind1 = intersect(z_behind,x_bc);
indices_behind = intersect(indices_behind1,y_bc);

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

[M_hcb,K_hcb,T_hcb,L_hcb] = CraigBamptonReduction(PrimalSub,200);

[V0s_hcb,om_hcb] = eigs(K_hcb,M_hcb, n_VMs, 'SM');
[f0_hcb,ind] = sort(sqrt(diag(om_hcb))/2/pi);
V0s_hcb = V0s_hcb(:,ind);

V0s_hcb_corr = converter_reducted_vector(PrimalSub,T_hcb,L_hcb,V0s_hcb(:,2));


V0_hcb = L_to_global(PrimalSub,V0s_hcb_corr);

corr_gmsh_indices = corr_gmsh_indices(PrimalSub,Mesh_ref);
V0_hcb = reindex_vector(corr_gmsh_indices,V0_hcb);

figure
hold on
nodalDef_hcb = reshape(V0_hcb,3,[]).';
PlotFieldonDeformedMesh(nodes_ref, elements_ref, nodalDef_hcb, 'factor', 1)
colormap jet
title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0_hcb(mod),3) ...
    ' Hz with HCB reduction'])

%% FRFs

S = 1e6;

%Applying a force on the substructured system
Fs = cell(1,PrimalSub.nSubs);

loc = [0.15 2.5 0.5];
locDOFs = PrimalSub.Substructures(2).Mesh.get_DOF_from_location(loc);
f=zeros(myMesh2.nDOFs,1);
f(locDOFs(1)) = S;
Fs{2} = f;

Fs{1} = zeros(myMesh1.nDOFs,1);
Fs{3} = zeros(myMesh3.nDOFs,1);
Fs{4} = zeros(myMesh4.nDOFs,1);
Fs{5} = zeros(myMesh5.nDOFs,1);

%Applying the force on the global system
locDOFs = Mesh_ref.get_DOF_from_location(loc);
F_ref=zeros(Mesh_ref.nDOFs,1);
F_ref(locDOFs(1)) = S;

Fc_ref = Assembly_ref.constrain_vector(F_ref);

om = sort(sqrt(diag(om)));
om_ref = sort(sqrt(diag(om_ref)));
om_hcb = sort(sqrt(diag(om_hcb)));

ksis = 0.1*ones(5,1);
A = [ones(5,1)./om/2 om.*ones(5,1)/2];
least_squares = (A'*A)\A'*ksis;
PrimalSub.DATA.C = least_squares(1)*PrimalSub.DATA.M + least_squares(2)*PrimalSub.DATA.K;
PrimalSub.DATA.Cc = PrimalSub.constrain_matrix(PrimalSub.DATA.C);

A_ref = [ones(5,1)./om_ref/2 om_ref.*ones(5,1)/2];
least_squares_ref = (A_ref'*A_ref)\A_ref'*ksis;
C_ref = least_squares(1)*M_ref + least_squares(2)*K_ref;
Cc_ref = Assembly_ref.constrain_matrix(C_ref);

A_hcb = [ones(5,1)./om_hcb/2 om_hcb.*ones(5,1)/2];
least_squares_hcb = (A_hcb'*A_hcb)\A_hcb'*ksis;
C_hcb = least_squares_hcb(1)*M_hcb + least_squares_hcb(2)*K_hcb;

fn = 100;
fmax = 200;


loc = [-2.5 0.15 0.5];
DOFs_s = PrimalSub.Substructures(1).Mesh.get_DOF_from_location(loc);
DOFs_s = DOFs_s(2);
DOFs_ref = Assembly_ref.Mesh.get_DOF_from_location(loc);
DOFs_ref = DOFs_ref(2);


[FRF] = frf_substructuring(PrimalSub,Fs,DOFs_s,fmax,fn);


%Calculating the FRF og the global model
FRF_ref = zeros(Mesh_ref.nDOFs,fn);
freqs = 0:fmax/(fn-1):fmax;
omegas = 2*pi*freqs;
for i=1:fn
    G = -omegas(i)^2*Mc_ref +1i*omegas(i)*Cc_ref +Kc_ref;
    FRF_ref(:,i) = Assembly_ref.unconstrain_vector(G\Fc_ref);
end


for xi=DOFs_ref
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

[FRF_hcb] = frf_hcb(PrimalSub,M_hcb,C_hcb,K_hcb,T_hcb,L_hcb,Fs,DOFs_s,fmax,fn);



%% Differences

diff_freq = abs(f0_ref-f0_hcb)./f0_ref;

normV0_ref = V0_ref(:,2)/norm(V0_ref(:,2));

normV0_hcb = V0_hcb/norm(V0_hcb);

hcb_angle = acos(normV0_hcb'*normV0_ref)*180/pi;

FRF_dof = abs(FRF(DOFs_s,:));
FRF_dof_ref = abs(FRF_ref(DOFs_ref,:));
FRF_dof_hcb = abs(FRF_hcb(DOFs_s,:));

normFRF = FRF_dof/norm(FRF_dof);
normFRF_ref = FRF_dof_ref/norm(FRF_dof_ref);
normFRF_hcb = FRF_dof_hcb/norm(FRF_dof_hcb);

FRF_sub_angle = acos(normFRF*normFRF_ref')*180/pi;
FRF_hcb_angle = acos(normFRF_hcb*normFRF_ref')*180/pi;







