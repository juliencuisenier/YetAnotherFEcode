function [Submeshes,Subelem] = Submeshing(Mesh,nSubs)

%Give as an output a cell containing Mesh structs of the substructures

Submeshes = {};

Subelem = {};

for iSub=1:nSubs
    
    iElements = [];

    for iElm=1:Mesh.nbElm
        
        if (Mesh.ELE_TAGS(iElm,1)==(20+iSub))||(Mesh.ELE_TAGS(iElm,1)==(20+iSub+2))...
                ||(Mesh.ELE_TAGS(iElm,1)==(20+iSub+4)) %Seeking if the j-th 
            %element belongs to the i-th substructure
            
            iElements = [iElements ; nonzeros(Mesh.ELE_NODES(iElm,:))];
            
        end

        iElements = unique(iElements); %Getting away of redunctant indices
        
        iNodes = Mesh.POS(iElements,:);
    
    end
    
    Submeshes{iSub} = iNodes;

end

end
