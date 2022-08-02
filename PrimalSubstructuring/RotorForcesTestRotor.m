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

%% SUBSTRUCTURING                                                   

% PrimalSubstructuring is a class with the basics of the substructuring and
% uses the primal assembly model of resolution 

globalIndices = get_globalIndices(Submeshes);
PrimalSub = PrimalSubstructuring([Assembly1 Assembly2...
    Assembly3 Assembly4 Assembly5],globalIndices,[]);

% Damping _________________________________________________________________
alfa = 67.2285;
beta = 1.3801e-4;
C = alfa*PrimalSub.DATA.M + beta*PrimalSub.DATA.K; % Rayleigh damping
PrimalSub.DATA.C = C;
PrimalSub.DATA.Cc = PrimalSub.constrain_matrix(C);


%% Rotor forces
Omega = [0 0 0];


PrimalSub.global_spinning_matrices(Omega);

n_VMs = 5;


[V0,om] = eigs(PrimalSub.DATA.Kc,PrimalSub.DATA.Mc, n_VMs, 'SM');
[f0,ind] = sort(sqrt(diag(om))/2/pi);
V0 = V0(:,ind);

V0  = PrimalSub.unconstrain_vector(V0);
V0s = L_to_local(PrimalSub,V0);


for mod=1:5
    figure
    hold on
    for jSub=1:PrimalSub.nSubs
        
        nodalDef = reshape(abs(V0s{jSub}(:,mod)),3,[]).';
        jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
        jElements = PrimalSub.Elements{jSub};
        PlotFieldonDeformedMesh(jMesh, jElements, nodalDef, 'factor', 10)
        
    end
    colormap jet
    title(['\Phi_' num2str(mod) ' - Frequency = ' num2str(f0(mod),3) ...
        ' Hz with substructuring'])
end