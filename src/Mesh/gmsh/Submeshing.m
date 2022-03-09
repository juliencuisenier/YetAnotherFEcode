function [Submeshes] = Submeshing(Mesh,nSubs)

%Give as an output a cell containing Mesh structs of the substructures such
%that : Submeshes{iSub,1} = array containing the coordinates of the nodes
%       Submeshes{iSub,2} = elements array
%       Submeshes{iSub,3} = array containing the original indices of the
%                           global mesh

Submeshes = {};


for iSub=1:nSubs
    
    
    iElements = [];

    for iElm=1:Mesh.nbHexas
        
        if Mesh.HEXAS(iElm,9)==24+iSub %Seeking if the j-th 
            %element belongs to the i-th substructure
            
            %Taking the element
            iElements = [iElements ; Mesh.ELE_NODES(iElm+1048,:)];
           
        end
        
    end
    
    iIndices = unique(iElements);
    
    iNodes = Mesh.POS(iIndices,:);
    
    Submeshes{iSub,1} = iNodes;
    
    Submeshes{iSub,2} = iElements;
    
    Submeshes{iSub,3} = iIndices;
    
end

end
