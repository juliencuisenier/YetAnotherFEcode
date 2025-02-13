function [V0,f0] = vibration_modes(PrimalSub,n_VMs,Mods,scaleFactor)
% Calculate and plot the vibration modes and natural frequencies of the
% substructured system
%- n_VMs is the number of VM that we want to considerate
% ex : n_VMs = 5 => we calculate the first 5 VMs
%- Mods is an array where there are the modes that we want to plot
%- scaleFactor is the scale factor for the plots
            
            [V0,om] = eigs(PrimalSub.DATA.Kc,PrimalSub.DATA.Mc, n_VMs, 'SM');
            [f0,ind] = sort(sqrt(diag(om))/2/pi);
            V0 = V0(:,ind);
            
            V0  = PrimalSub.unconstrain_vector(V0);
            
            for iMod=Mods
                
                figure
                hold on
                
                for jSub=1:PrimalSub.nSubs
                    
                    nodalDef = reshape(V0{jSub}(:,iMod),3,[]).';
                    jMesh = PrimalSub.Substructures(jSub).Mesh.nodes;
                    jElements = PrimalSub.Elements{jSub};
                    PlotFieldonDeformedMesh(jMesh, jElements, nodalDef, 'factor', scaleFactor)
         
                end
                colormap jet
                title(['\Phi_' num2str(iMod) ' - Frequency = ' num2str(f0(iMod),3) ...
                ' Hz with substructuring'])
            end
        end
