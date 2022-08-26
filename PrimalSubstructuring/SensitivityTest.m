% Initialize Matlab
clear
close all
clc
format shortg


%% SETTINGS                                                         

% if USEJULIA == 1, tensors are computed with Julia, otherwise with Matlab.
% If you set USEJULIA to 1 make sure you have Julia and Mex.jl installed
% (see https://github.com/byuflowlab/Mex.jl)
USEJULIA = 0;

% Sensitivity Analysis settings ___________________________________________
SENS = 1;	% perform sensitivity analysis (1) or used stored results (0).
filename_sens = 'sens_2VM_H3';	% file with sens analysis results stored
method_sens = 'normal';         % available options are 'normal'/'omegaconst'

% number of vibration modes used in modal reduction (MAKE SURE it's the
% same as in the loaded data, if loading sensitivities - see above)
n_VMs = 5;

% do you want to save sensitivity analysis results?
save_sens = 0;
save_as = 'sens_2VM_H3.mat';	% name of the file you save results

% DpROM settings __________________________________________________________
% parametric formulation for defects
FORMULATION = 'N1'; % N1/N1t/N0
VOLUME = 1;         % integration over defected (1) or nominal volume (0)

% defect amplitudes
xi1 = 0.2;     % defect amplitude - top blade angle (x)
xi2 = 0.3;     % defect amplitude - top blade angle (z)
xi3 = 0.3;     % defect amplitude - bottom blade angle (y)
xi4 = 0.4;     % defect amplitude - left blade angle (z)

% Common settings for HB __________________________________________________
imod = 1;           % eigenfreq to study
H = 3;             	% harmonic order
startfactor = 1.2;  % start frequency = f0(imod) * startfactor
endfactor = 0.8;    %   end frequency = f0(imod) * endfactor
ds = 20;           	% Path continuation step size
exc_lev = 1;     % excitation level


%% Prepare Model (full nominal, full defected)                      

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = 0.3; % thickness of a blade
L = 2; % length of a blade
Omega = [0 0 147]; % rotation speed in rad/s

% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = false;	% set "true" for plane_stress
% Element
myElementConstructor = @()Hex8Element(myMaterial);

% MESH_____________________________________________________________________

TotalMesh = load_gmsh2("RotorSketch.msh", -1);
[nodes,elements,nset] = extract_gmsh(TotalMesh);

% nominal mesh
NominalMesh = Mesh(nodes);
NominalMesh.create_elements_table(elements,myElementConstructor);
NominalMesh.set_essential_boundary_condition(nset{3},1:3,0)
NominalMesh.set_essential_boundary_condition(nset{6},1:3,0)

% DEFECT SHAPES ***********************************************************
% (1) angle defect on the top blade (displacement along x axis)
nodes_r = zeros(NominalMesh.nNodes,1);
nodes_r(nodes(:,2)>0.5) = 1;
th = 5; %angle in �
ud = (nodes(:,2))*tan(th*pi/180).*nodes_r;
angle_defect_top_x = zeros(numel(nodes),1);
angle_defect_top_x(1:3:end) = ud;

% (2) angle defect on the top blade (displacement along z axis)
nodes_r = zeros(NominalMesh.nNodes,1);
nodes_r(nodes(:,2)>0.5) = 1;
th2 = 5; %angle in �
ud = (nodes(:,2))*tan(th2*pi/180).*nodes_r;
angle_defect_top_z = zeros(numel(nodes),1);
angle_defect_top_z(3:3:end) = ud;

% (3) angle defect on the right blade (displacement along y axis)
nodes_r = zeros(NominalMesh.nNodes,1);
nodes_r(nodes(:,1)>0.5) = 1;
th3 = 5; %angle in �
ud = (nodes(:,1))*tan(th3*pi/180).*nodes_r;
angle_defect_right_y = zeros(numel(nodes),1);
angle_defect_right_y(1:3:end) = ud;

% (4) angle defect on the left blade (displacement along z-axis)
nodes_r = zeros(NominalMesh.nNodes,1);
nodes_r(nodes(:,1)<-0.5) = 1;
th4 = 5; %angle in �
ud = (nodes(:,1))*tan(th4*pi/180).*nodes_r;
angle_defect_left_z = zeros(numel(nodes),1);
angle_defect_left_z(3:3:end) = ud;

% *************************************************************************


% defected mesh
 U = [angle_defect_top_x,angle_defect_top_z,...
     angle_defect_right_y, angle_defect_left_z]; % defect basis
xi = [xi1 xi2 xi3 xi4]'; % parameter vector
m = length(xi); % number of parameters

% update defected mesh nodes
d = U*xi; % displacement fields introduced by defects
dd = [d(1:3:end) d(2:3:end) d(3:3:end)] ;
nodes_defected = nodes + dd; % nominal + d ---> defected 
DefectedMesh = Mesh(nodes_defected);
DefectedMesh.create_elements_table(elements,myElementConstructor);
DefectedMesh.set_essential_boundary_condition(nset{3},1:3,0)
DefectedMesh.set_essential_boundary_condition(nset{6},1:3,0)
    fprintf(' Roughness defect (x direction):         %.2f * sin(100*pi*x/L) \n', xi1)
    fprintf(' Roughness defect (z direction):         %.2f * sin(100*pi*z/L) \n', xi2)
    fprintf(' Top blade angle defect (x direction):         %.1f째 (xi = %.1f) \n', th*xi1, xi1)
    fprintf(' Top blade angle defect (z direction):         %.1f째 (xi = %.1f) \n', th2*xi2, xi2)
    fprintf(' Right blade angle defect (y direction):         %.1f째 (xi = %.1f) \n', th3*xi3, xi3)
    fprintf(' Left blade angle defect (z direction):         %.1f째 (xi = %.1f) \n', th4*xi4, xi4)
figure('units','normalized','position',[.2 .3 .6 .4])
v1 = reshape(U*xi, 3, []).';
S = 1;%2 * Ly / max(abs(v1(:)));
hf=PlotFieldonDeformedMesh(nodes, elements, v1, 'factor', S);
% title(sprintf('Defect, \\xi=[ %.2f, %.2f, %.lf, %.1f ], S=%.1f\\times',...
%     xi1, xi2, xi3, xi4, S))
%axis normal; grid on; box on; set(hf{1},'FaceAlpha',.7); drawnow



% ASSEMBLY ________________________________________________________________
% nominal
NominalAssembly = Assembly(NominalMesh);
Mn = NominalAssembly.mass_matrix();
nNodes = size(nodes, 1);
u0 = zeros( NominalMesh.nDOFs, 1);
[Kn,~] = NominalAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    NominalAssembly.DATA.K = Kn;
    NominalAssembly.DATA.M = Mn;
    
% Damping _________________________________________________________________
alfa = 67.2285;
beta = 1.3801e-4;
D = alfa*Mn + beta*Kn; % Rayleigh damping
NominalAssembly.DATA.D = D;

NominalG = NominalAssembly.coriolis_matrix(Omega);
NominalKsp = NominalAssembly.spin_softening_matrix(Omega);

NominalAssembly.DATA.D = NominalAssembly.DATA.D + NominalG;
NominalAssembly.DATA.K = NominalAssembly.DATA.K + NominalKsp;
NominalAssembly.DATA.Dc = NominalAssembly.constrain_matrix(NominalAssembly.DATA.D);
NominalAssembly.DATA.Kc = NominalAssembly.constrain_matrix(NominalAssembly.DATA.K);

% defected
DefectedAssembly = Assembly(DefectedMesh);
Md = DefectedAssembly.mass_matrix();
u0 = zeros( DefectedMesh.nDOFs, 1);
[Kd,~] = DefectedAssembly.tangent_stiffness_and_force(u0);
    % store matrices
    DefectedAssembly.DATA.K = Kd;
    DefectedAssembly.DATA.M = Md;
    
    % Damping _________________________________________________________________
DefectedD = alfa*Md + beta*Kd; % Rayleigh damping
DefectedAssembly.DATA.D = DefectedD;

DefectedG = DefectedAssembly.coriolis_matrix(Omega);
DefectedKsp = DefectedAssembly.spin_softening_matrix(Omega);

DefectedAssembly.DATA.D = DefectedAssembly.DATA.D + DefectedG;
DefectedAssembly.DATA.K = DefectedAssembly.DATA.K + DefectedKsp;
DefectedAssembly.DATA.Dc = DefectedAssembly.constrain_matrix(DefectedAssembly.DATA.D);
DefectedAssembly.DATA.Kc = DefectedAssembly.constrain_matrix(DefectedAssembly.DATA.K);

% External force __________________________________________________________
% blades are forced clockwise
S = 1e6;

Fext = zeros(NominalMesh.nDOFs,1);

% Force on left blade
loc1 = [-2.5 -0.15 0.5];
DOFs = NominalMesh.get_DOF_from_location(loc1);
Fext(DOFs(2)) = S;

% Force on top blade
loc2 = [-0.15 2.5 0.5];
DOFs = NominalMesh.get_DOF_from_location(loc2);
Fext(DOFs(1)) = S;

% Force on right blade
loc3 = [2.5 0.15 0.5];
DOFs = NominalMesh.get_DOF_from_location(loc3);
Fext(DOFs(2)) = -S;

% Force on bottom blade
loc4 = [0.15 -2.5 0.5];
DOFs = NominalMesh.get_DOF_from_location(loc4);
Fext(DOFs(1)) = -S;


%% Eigenmodes, number of vibration modes in ROM                     

% Eigenvalue problem_______________________________________________________
% Vibration Modes (VM): nominal
Mnc = NominalAssembly.constrain_matrix(Mn);
[VMn,om] = eigs(NominalAssembly.DATA.Kc, Mnc, n_VMs, 'SM');
[f0n,ind] = sort(sqrt(diag(om))/2/pi);
VMn = VMn(:,ind);
for ii = 1:n_VMs
    VMn(:,ii) = VMn(:,ii)/norm(VMn(:,ii));
end
VMn = NominalAssembly.unconstrain_vector(VMn);

% Vibration Modes (VM): defected
Mdc = DefectedAssembly.constrain_matrix(Md);
[VMd,om] = eigs(DefectedAssembly.DATA.Kc, Mdc, n_VMs, 'SM');
[f0d,ind] = sort(sqrt(diag(om))/2/pi);
VMd = VMd(:,ind);
for ii = 1:n_VMs
    VMd(:,ii) = VMd(:,ii)/norm(VMd(:,ii));
end
VMd = DefectedAssembly.unconstrain_vector(VMd);

angles_VMs = zeros(1,n_VMs);

for ii=1:n_VMs
    angles_VMs(ii) = acos(VMd(:,ii)'*VMn(:,ii))*180/pi;
end



%% FRFs 

Fextc = NominalAssembly.constrain_vector(Fext);

fmin = 0 ; %minimal frequency
fmax = 200 ; % maximal frequency
fn = 50 ; % number of point of the FRF's curve

NominalFRF = zeros(NominalMesh.nDOFs,fn);
DefectedFRF = zeros(NominalMesh.nDOFs,fn);

freqs = fmin:(fmax-fmin)/(fn-1):fmax;
omegas = 2*pi*freqs;
for i=1:fn
    G = -omegas(i)^2*Mnc +1i*omegas(i)*NominalAssembly.DATA.Dc + NominalAssembly.DATA.Kc;
    NominalFRF(:,i) = NominalAssembly.unconstrain_vector(G\Fextc);
    G = -omegas(i)^2*Mdc +1i*omegas(i)*DefectedAssembly.DATA.Dc + DefectedAssembly.DATA.Kc;
    DefectedFRF(:,i) = DefectedAssembly.unconstrain_vector(G\Fextc);
end

loc = [0.15 1.5 0.5];
DOFs = NominalMesh.get_DOF_from_location(loc);


for x_i=DOFs
   figure
   hold on
   subplot(2,1,1);
   plot(freqs,abs(NominalFRF(x_i,:)),freqs,abs(DefectedFRF(x_i,:)));
   title(strcat("Module of z-DOF's FRF"));
   legend('Nominal model','Defected model');
   xlabel('Frequency');
   set(gca,'yscale','log');
   subplot(2,1,2);
   plot(freqs,unwrap(angle(NominalFRF(x_i,:))),freqs,unwrap(angle(DefectedFRF(x_i,:))))
   title(strcat("Phase of the DOF z-DOF's FRF [rad]"));
   xlabel('Frequency');
   legend('Nominal model','Defected model');
end

% %% Modal Derivatives & Defect Sensitivities                         
% 
% % nominal
% [MDn, MDname] = modal_derivatives(NominalAssembly, elements, VMn);
% % defected
% MDd = modal_derivatives(DefectedAssembly, elements, VMd);
% 
% % defect sensitivities
% [DS, names] = defect_sensitivities(NominalAssembly, elements, VMn, U, ...
%     FORMULATION);
% 
% %% Generate ROM tensors (ROM-n DpROM ROM-d)                         
% % define reduced order basis
% Vn = [VMn MDn];     % reduced order basis (ROM-n)
% V  = [VMn MDn DS]; 	% reduced order basis (DpROM)
% Vd = [VMd MDd];   	% reduced order basis (ROM-d)
% 
% % orthonormalize reduction basis
% Vn = orth(Vn);	% ROM-n
% V  = orth(V);	% DpROM
% Vd = orth(Vd);	% ROM-d
% 
% 
% % standard reduced order model (no defects in the mesh)
% tensors_ROMn = reduced_tensors_ROM(NominalAssembly, elements, Vn, USEJULIA);
% 
% % standard reduced order model (defects in the mesh)
% tensors_ROMd = reduced_tensors_ROM(DefectedAssembly, elements, Vd, USEJULIA);
% tensors_ROMd.xi = xi; % save for which xi ROMd is computed
% 
% % parametric formulation for defects
% tensors_DpROM = reduced_tensors_DpROM(NominalAssembly, elements, ...
%     V, U, FORMULATION, VOLUME, USEJULIA); %compute tensors
% 
% % evaluate the defected tensors at xi
% [Q2,Mxi] = DefectedTensors(tensors_DpROM, xi);
% 
% Mnr = Vn'*Mn*Vn; 	% reduced mass matrix (ROM-n)
% Mdr = Vd'*Md*Vd; 	% reduced mass matrix (ROM-d)
% if VOLUME == 0
%     Mr  = V' *Mn*V; 	% reduced mass matrix (DpROM)
% elseif VOLUME == 1
%     Mr  = V' *Md*V; 	% reduced mass matrix (DpROM)
% end
% 
% % compute eigenfrequencies of reduced models
% f0_ROMn = sort(sqrt( eigs(tensors_ROMn.Q2, Mnr, n_VMs, 'SM') )/2/pi);
% f0_ROMd = sort(sqrt( eigs(tensors_ROMd.Q2, Mdr, n_VMs, 'SM') )/2/pi); 
% f0_DpROM = sort(sqrt( eigs(Q2, Mr, n_VMs, 'SM') )/2/pi);
% f0_DpROM_Mxi = sort(sqrt( eigs(Q2, Mxi, n_VMs, 'SM') )/2/pi);


%% %%%%%%%%%%%%%%%%% Substructuring %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparing the substructured nominal model

% SUBMESHING:______________________________________________________________

TotalMesh = load_gmsh2("RotorSketch.msh", -1); %reading the .msh file
Submeshes = Submeshing(TotalMesh);


%Mesh 1 (Left blade)_______________________________________________________

NominalMesh1 = Mesh(Submeshes{1,1});
NominalMesh1.create_elements_table(Submeshes{1,2}, myElementConstructor);

%Mesh 2 (Top blade)________________________________________________________

NominalMesh2 = Mesh(Submeshes{2,1});
NominalMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);


%Mesh 3 (Right blade)______________________________________________________

NominalMesh3 = Mesh(Submeshes{3,1});
NominalMesh3.create_elements_table(Submeshes{3,2}, myElementConstructor);



%Mesh 4 (Bottom blade)_____________________________________________________

NominalMesh4 = Mesh(Submeshes{4,1});
NominalMesh4.create_elements_table(Submeshes{4,2}, myElementConstructor);


%Mesh 5 (Heart)____________________________________________________________

NominalMesh5 = Mesh(Submeshes{5,1});
NominalMesh5.create_elements_table(Submeshes{5,2}, myElementConstructor);

NominalMesh5.set_essential_boundary_condition(Submeshes{5,3}{3},1:3,0)
NominalMesh5.set_essential_boundary_condition(Submeshes{5,3}{6},1:3,0)


% ASSEMBLY (nominal)_____________________________________________________
NominalAssembly1 = Assembly(NominalMesh1);
NominalAssembly2 = Assembly(NominalMesh2);
NominalAssembly3 = Assembly(NominalMesh3);
NominalAssembly4 = Assembly(NominalMesh4);
NominalAssembly5 = Assembly(NominalMesh5);
                                                 
globalIndices = get_globalIndices(Submeshes);
NominalPrimalSub = PrimalSubstructuring([NominalAssembly1 NominalAssembly2...
    NominalAssembly3 NominalAssembly4 NominalAssembly5],globalIndices,[]);

% Damping _________________________________________________________________
% The Rayleigh damping coefficients are the same as the full global model
NominalPrimalSub.DATA.C = alfa*NominalPrimalSub.DATA.M + beta*NominalPrimalSub.DATA.K;
NominalPrimalSub.DATA.Cc = NominalPrimalSub.constrain_matrix(NominalPrimalSub.DATA.C);

% Rotor forces_____________________________________________________________

NominalPrimalSub.global_spinning_matrices(Omega);

%% DEFECT SHAPES ***********************************************************

%(i) Top blade

nodes = NominalMesh2.nodes;

% (1) roughness on the top blade (along x)
nodes_r = zeros(NominalMesh2.nNodes,1);
nodes_r(nodes(:,2)>0.55) = 1; % selecting the nodes above y=0.5
vd1 = sin(100*pi/L * nodes(:,1).*nodes_r);	
roughness_defect_x_sub = zeros(numel(nodes),1);
roughness_defect_x_sub(1:3:end) = vd1; % vectorized defect-field on the substructure
roughness_defect_x = NominalPrimalSub.L{2}'*roughness_defect_x_sub; 


% (2) roughness on the top blade (along z)
nodes_r = zeros(NominalMesh2.nNodes,1);
nodes_r(nodes(:,2)>0.55) = 1; % selecting the nodes above y=0.5
vd1 = sin(100*pi/L * nodes(:,3).*nodes_r);	
roughness_defect_z_sub = zeros(numel(nodes),1);
roughness_defect_z_sub(3:3:end) = vd1; % vectorized defect-field
roughness_defect_z = NominalPrimalSub.L{2}'*roughness_defect_z_sub;

% (3) angle defect on the top blade (displacement along x axis)
nodes_r = zeros(NominalMesh2.nNodes,1);
nodes_r(nodes(:,2)>0.5) = 1;
th = 5; %angle in �
ud = (nodes(:,2))*tan(th*pi/180).*nodes_r;
angle_defect_top_x_sub = zeros(numel(nodes),1);
angle_defect_top_x_sub(1:3:end) = ud;
angle_defect_top_x = NominalPrimalSub.L{2}'*angle_defect_top_x_sub;

% (4) angle defect on the top blade (displacement along z axis)
nodes_r = zeros(NominalMesh2.nNodes,1);
nodes_r(nodes(:,2)>0.5) = 1;
th2 = 5; %angle in �
ud = (nodes(:,2))*tan(th2*pi/180).*nodes_r;
angle_defect_top_z_sub = zeros(numel(nodes),1);
angle_defect_top_z_sub(3:3:end) = ud;
angle_defect_top_z = NominalPrimalSub.L{2}'*angle_defect_top_z_sub;

%(ii) Right blade

nodes = NominalMesh3.nodes;

% (5) angle defect on the right blade (displacement along y axis)
nodes_r = zeros(NominalMesh3.nNodes,1);
nodes_r(nodes(:,1)>0.5) = 1;
th3 = 5; %angle in �
ud = (nodes(:,1))*tan(th3*pi/180).*nodes_r;
angle_defect_right_y_sub = zeros(numel(nodes),1);
angle_defect_right_y_sub(1:3:end) = ud;
angle_defect_right_y = NominalPrimalSub.L{3}'*angle_defect_right_y_sub;

% (iii) Left blade

nodes = NominalMesh1.nodes;

% (6) angle defect on the left blade (displacement along z axis)

nodes_r = zeros(NominalMesh3.nNodes,1);
nodes_r(nodes(:,1)>0.5) = 1;
th4 = 5; %angle in �
ud = (nodes(:,1))*tan(th4*pi/180).*nodes_r;
angle_defect_left_z_sub = zeros(numel(nodes),1);
angle_defect_left_z_sub(3:3:end) = ud;
angle_defect_left_z = NominalPrimalSub.L{1}'*angle_defect_left_z_sub;
% *************************************************************************


% defected mesh
 U = [roughness_defect_x, roughness_defect_z, angle_defect_top_x,...
     angle_defect_top_z, angle_defect_right_y,angle_defect_left_z]; % defect basis

% update defected mesh nodes
d = U*xi; % displacement fields introduced by defects

%% Preparing the substructured defected model

% (i) Left blade

d_sub = NominalPrimalSub.L{1}*d;
dd = [d_sub(1:3:end) d_sub(2:3:end) d_sub(3:3:end)];
nodes_defected = NominalMesh1.nodes + dd;
DefectedMesh1 = Mesh(nodes_defected);
DefectedMesh1.create_elements_table(Submeshes{1,2}, myElementConstructor);

% (ii) Top blade

d_sub = NominalPrimalSub.L{2}*d;
dd = [d_sub(1:3:end) d_sub(2:3:end) d_sub(3:3:end)];
nodes_defected = NominalMesh2.nodes + dd;
DefectedMesh2 = Mesh(nodes_defected);
DefectedMesh2.create_elements_table(Submeshes{2,2}, myElementConstructor);

% (iii) Right blade

d_sub = NominalPrimalSub.L{3}*d;
dd = [d_sub(1:3:end) d_sub(2:3:end) d_sub(3:3:end)];
nodes_defected = NominalMesh3.nodes + dd;
DefectedMesh3 = Mesh(nodes_defected);
DefectedMesh3.create_elements_table(Submeshes{3,2}, myElementConstructor);

% (iv) Bottom blade

d_sub = NominalPrimalSub.L{4}*d;
dd = [d_sub(1:3:end) d_sub(2:3:end) d_sub(3:3:end)];
nodes_defected = NominalMesh4.nodes + dd;
DefectedMesh4 = Mesh(nodes_defected);
DefectedMesh4.create_elements_table(Submeshes{4,2}, myElementConstructor);

% (v) Heart

d_sub = NominalPrimalSub.L{5}*d;
dd = [d_sub(1:3:end) d_sub(2:3:end) d_sub(3:3:end)];
nodes_defected = NominalMesh5.nodes + dd;
DefectedMesh5 = Mesh(nodes_defected);
DefectedMesh5.create_elements_table(Submeshes{5,2}, myElementConstructor);

DefectedMesh5.set_essential_boundary_condition(Submeshes{5,3}{3},1:3,0)
DefectedMesh5.set_essential_boundary_condition(Submeshes{5,3}{6},1:3,0)

% Assembly (defected)

DefectedAssembly1 = Assembly(DefectedMesh1);
DefectedAssembly2 = Assembly(DefectedMesh2);
DefectedAssembly3 = Assembly(DefectedMesh3);
DefectedAssembly4 = Assembly(DefectedMesh4);
DefectedAssembly5 = Assembly(DefectedMesh5);

DefectedPrimalSub = PrimalSubstructuring([DefectedAssembly1 DefectedAssembly2...
    DefectedAssembly3 DefectedAssembly4 DefectedAssembly5],globalIndices,[]);

% Damping _________________________________________________________________
DefectedPrimalSub.DATA.C = alfa*DefectedPrimalSub.DATA.M + beta*DefectedPrimalSub.DATA.K;
DefectedPrimalSub.DATA.Cc = DefectedPrimalSub.constrain_matrix(DefectedPrimalSub.DATA.C);

% Spinning forces__________________________________________________________

DefectedPrimalSub.global_spinning_matrices(Omega);

% External force __________________________________________________________
% blades are forced clockwise

% Force on left blade
f1 = zeros(NominalMesh1.nDOFs,1);
% The location are the same as full model
DOFs = NominalMesh1.get_DOF_from_location(loc1);
f1(DOFs(2)) = S;

% Force on top blade
f2 = zeros(NominalMesh2.nDOFs,1);
DOFs = NominalMesh2.get_DOF_from_location(loc2);
f2(DOFs(1)) = S;

% Force on right blade
f3 = zeros(NominalMesh3.nDOFs,1);
DOFs = NominalMesh3.get_DOF_from_location(loc3);
f3(DOFs(2)) = -S;

% Force on bottom blade
f4 = zeros(NominalMesh4.nDOFs,1);
DOFs = NominalMesh4.get_DOF_from_location(loc4);
f4(DOFs(1)) = -S;

% Force on heart
f5 = zeros(NominalMesh5.nDOFs,1); % no forces applied on it

Fext = {f1,f2,f3,f4,f5};

%% Eigenmodes, number of vibration modes in ROM                     

% Eigenvalue problem_______________________________________________________
% Vibration Modes (VM): nominal
[VMn_sub,om] = eigs(NominalPrimalSub.DATA.Kc, NominalPrimalSub.DATA.Mc, n_VMs, 'SM');
[f0n_sub,ind] = sort(sqrt(diag(om))/2/pi);
VMn_sub = VMn_sub(:,ind);
for ii = 1:n_VMs
    VMn_sub(:,ii) = VMn_sub(:,ii)/max(sqrt(sum(VMn_sub(:,ii).^2,2)));
end
VMn_sub = NominalPrimalSub.unconstrain_vector(VMn_sub);

% Vibration Modes (VM): defected
[VMd_sub,om] = eigs(DefectedPrimalSub.DATA.Kc, DefectedPrimalSub.DATA.Mc, n_VMs, 'SM');
[f0d_sub,ind] = sort(sqrt(diag(om))/2/pi);
VMd_sub = VMd_sub(:,ind);
for ii = 1:n_VMs
    VMd_sub(:,ii) = VMd_sub(:,ii)/max(sqrt(sum(VMd_sub(:,ii).^2,2)));
end
VMd_sub = DefectedPrimalSub.unconstrain_vector(VMd_sub);

%% FRFs 
Fext_sub = zeros(NominalPrimalSub.nDOFglobal,1);

for iSub = 1:NominalPrimalSub.nSubs
    Fext_sub = Fext_sub + NominalPrimalSub.L{iSub}'*Fext{iSub};
end

Fextc_sub = NominalPrimalSub.constrain_vector(Fext_sub);

loc = [0.15 1.5 0.5];
DOFs = NominalMesh2.get_DOF_from_location(loc);
DOFs_sub = NominalPrimalSub.Us{2}(DOFs);

fmin = 0 ; %minimal frequency
fmax = 200 ; % maximal frequency
fn = 50 ; % number of point of the FRF's curve

% NominalFRF = zeros(NominalPrimalSub.nDOFglobal,fn);
% DefectedFRF = zeros(NominalPrimalSub.nDOFglobal,fn);
% 
% freqs = fmin:(fmax-fmin)/(fn-1):fmax;
% omegas = 2*pi*freqs;
% for i=1:fn
%     G = -omegas(i)^2*NominalPrimalSub.DATA.Mc +1i*omegas(i)*NominalPrimalSub.DATA.Cc + NominalPrimalSub.DATA.Kc;
%     NominalFRF(:,i) = NominalPrimalSub.unconstrain_vector(G\Fextc_sub);
%     G = -omegas(i)^2*DefectedPrimalSub.DATA.Mc +1i*omegas(i)*DefectedPrimalSub.DATA.Cc + DefectedPrimalSub.DATA.Kc;
%     DefectedFRF(:,i) = DefectedPrimalSub.unconstrain_vector(G\Fextc_sub);
% end
% 
% 
% for x_i=DOFs
%    figure
%    hold on
%    subplot(2,1,1);
%    plot(freqs,abs(NominalFRF(x_i,:)),freqs,abs(DefectedFRF(x_i,:)));
%    title(strcat('Module of the DOF ',num2str(x_i),"'s FRF "));
%    legend('Nominal model','Defected model');
%    xlabel('Frequency');
%    set(gca,'yscale','log');
%    subplot(2,1,2);
%    plot(freqs,unwrap(angle(NominalFRF(x_i,:))),freqs,unwrap(angle(DefectedFRF(x_i,:))))
%    title(strcat('Phase of the DOF ',num2str(x_i),"'s FRF [rad]"));
%    xlabel('Frequency');
%    legend('Nominal model','Defected model');
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Performing Craig-Bampton reduction

% On Nominal
[Mn_cb,Cn_cb,Kn_cb,Tn_cb,Ln_cb] = CraigBamptonReduction(NominalPrimalSub,Omega,200);


Cn_cb = alfa*Mn_cb + beta*Kn_cb + Cn_cb;


[VMn_cb,om] = eigs(Kn_cb, Mn_cb, n_VMs, 'SM');
[f0n_cb,ind] = sort(sqrt(diag(om))/2/pi);


% On defected
[Md_cb,Cd_cb,Kd_cb,Td_cb,Ld_cb] = CraigBamptonReduction(DefectedPrimalSub,Omega,200);


Cd_cb = alfa*Md_cb + beta*Kd_cb + Cd_cb;

[VMd_cb,om] = eigs(Kd_cb, Md_cb, n_VMs, 'SM');
[f0d_cb,ind] = sort(sqrt(diag(om))/2/pi);

angles_VMs_cb = zeros(1,n_VMs);

corrector = corr_gmsh_indices(DefectedPrimalSub,DefectedMesh);

for ii = 1:n_VMs
    VMd_cb_loc = converter_reducted_vector(DefectedPrimalSub,Td_cb,Ld_cb,VMd_cb(:,ii));
    VMd_cb_loc = reindex_vector(corrector,VMd_cb_loc);
    
    VMd_cb_loc = VMd_cb_loc/norm(VMd_cb_loc);
    
    angles_VMs_cb(ii) = acos(VMd(:,ii)'*VMd_cb_loc)*180/pi;
end

[DefectedFRF_cb] = frf_cb(DefectedPrimalSub,Md_cb,Cd_cb,Kd_cb,Td_cb,Ld_cb,Fext,DOFs,fmin,fmax,fn);
[NominalFRF_cb] = frf_cb(NominalPrimalSub,Mn_cb,Cn_cb,Kn_cb,Tn_cb,Ln_cb,Fext,DOFs,fmin,fmax,fn);

angles_FRF = zeros(1,fn);


for ii=1:fn
    DefectedFRF_cb(:,ii) = reindex_vector(corrector,abs(DefectedFRF_cb (:,ii)));
    
    normFRF_cb = DefectedFRF_cb(:,ii)/norm(DefectedFRF_cb(:,ii));
    normFRF = abs(DefectedFRF(:,ii))/norm(abs(DefectedFRF(:,ii)));
    
    angles_FRF(ii) = acos(normFRF'*normFRF_cb)*180/pi;
end

% for ii=1:fn
%     NominalFRF_cb(:,ii) = reindex_vector(corrector,abs(NominalFRF_cb (:,ii)));
%     
%     normFRF_cb = NominalFRF_cb(:,ii)/norm(NominalFRF_cb(:,ii));
%     normFRF = abs(NominalFRF(:,ii))/norm(abs(NominalFRF(:,ii)));
%     
%     angles_FRF(ii) = acos(normFRF'*normFRF_cb)*180/pi;
% end


mean(angles_FRF)
