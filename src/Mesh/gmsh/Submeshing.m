function [Submeshes] = Submeshing(Mesh,nSubs)

%Give as an output a cell containing Mesh structs of the substructures

Submeshes = {};


for iSub=1:nSubs
    
    
    iElements = [];

    for iElm=1:Mesh.nbElm
        
        if (Mesh.ELE_TAGS(iElm,1)==(20+iSub))||(Mesh.ELE_TAGS(iElm,1)==(20+iSub+2))...
                ||(Mesh.ELE_TAGS(iElm,1)==(20+iSub+4)) %Seeking if the j-th 
            %element belongs to the i-th substructure
            
            %Taking the element
            iElements = [iElements ; Mesh.ELE_NODES(iElm,:)];
           
            
        end
        
        %Taking the nodes indices which belong to the iSub substructure
        iIndices = unique(nonzeros(iElements));
        
        iNodes = Mesh.POS(iIndices,:);
    
    end
    
    Submeshes{iSub,1} = iNodes;
    
    Submeshes{iSub,2} = iElements;
    
    Submeshes{iSub,3} = [iIndices];

end

end
