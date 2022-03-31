function [V0s,V0s_defect,f0] = defect_modes(PrimalSub,F,n_VMs,Mods)
%DEFECT_MODES gives the result and the plot of the original PrimalSub modes
%and the same for a defected PrimalSub, defect caused by the forces 'cell F 

%   PrimalSub : PrimalSubstructure class
%   F : cell of forces applied on each substructure
%   n_VMs : number of the first vibration mode that we want
%   Mods : Modes that we want to plot

[V0s,f0] = PrimalSub.vibration_mode(n_VMs,Mods);

u = PrimalSub.static_resolution(F);

us = L_to_local(PrimalSub,u);

PrimalSubDefect = PrimalSubstructuring(PrimalSub.Substructures,...
    PrimalSub.globalIndices,us);

V0s_defect = PrimalSubDefect.vibration_mode(n_VMs,Mods);

% [M_defect,K_defect] = PrimalSub.global_mass_stiffness(us);
% 
% Mc_defect = PrimalSub.constrain_matrix(M_defect);
% Kc_defect = PrimalSub.constrain_matrix(K_defect);
% 
% [V0,om] = eigs(Mc_defect,Kc_defect, n_VMs, 'SM');
% [f0,ind] = sort(sqrt(diag(om))/2/pi);
% V0 = V0(:,ind);
% 
% V0  = PrimalSub.unconstrain_vector(V0);
% 
% V0s_defect = L_to_local(PrimalSub,V0);
% 
% for iMod=Mods
%     
%     figure
%     hold on
%     
%     for jSub=1:PrimalSub.nSubs
%         
%         nodalDef = reshape(V0s_defect{jSub}(:,iMod),3,[]).';
%         jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
%         jElements = PrimalSub.Elements{jSub};
%         PlotFieldonDeformedMesh(jMesh, jElements, nodalDef, 'factor', 1)
%         
%     end
%     colormap jet
%     title(['\Phi_' num2str(iMod) ' - Frequency = ' num2str(f0(iMod),3) ...
%         ' Hz with defected substructuring'])
% end


end

