function [Submeshes] = Submeshing(Mesh)

%Give as an output a cell containing Mesh structs of the substructures such
%that : Submeshes{iSub,1} = array containing the coordinates of the nodes
%       Submeshes{iSub,2} = elements array
%       Submeshes{iSub,3} = indices of the boundaries
%       Submeshes{iSub,4} = array containing the original indices of the
%                           global mesh



nSubs = max(Mesh.ELE_TAGS(:,1))-min(Mesh.ELE_TAGS(:,1))+1;

Submeshes = cell(nSubs,4);

min_index = min(Mesh.ELE_TAGS(:,1)) - 1;

for iSub=1:nSubs
    
    %Find the indices of elements belonging to the i-th substructure
    iElementsIndices = Mesh.ELE_TAGS(:,1)==min_index+iSub;
    
    %Taking the elements
    iElements = Mesh.ELE_NODES(iElementsIndices,:);
    
    %Taking the indices of the nodes of the elements
    iIndices = unique(iElements);
    
    %Reindex the Elements array
    [nbElm, nbNodes_by_Elm] = size(iElements);  
    reindexed_elements = zeros(nbElm, nbNodes_by_Elm);
  
    for iElm=1:nbElm 
        for jNode=1:nbNodes_by_Elm
            new_indice = find(iIndices == iElements(iElm,jNode));
            reindexed_elements(iElm,jNode) = new_indice;      
        end  
    end
    iElements = reindexed_elements;
    
    %Taking the coordinates of the nodes which belong to the i-th
    %substructure
    iNodes = Mesh.POS(iIndices,:);
    
    %Finding the min and max for x,y and z
    iNset = cell(1,6);
    
    for j =1:3
        iNset{j} = find(iNodes(:,j)== min(iNodes(:,j)));
        iNset{3+j} = find(iNodes(:,j)== max(iNodes(:,j)));     
    end
    
    
    %Filling the Submeshes cell
    Submeshes{iSub,1} = iNodes;
    
    Submeshes{iSub,2} = iElements;
    
    Submeshes{iSub,3} = iNset;
    
    Submeshes{iSub,4} = iIndices;
    
end

end
