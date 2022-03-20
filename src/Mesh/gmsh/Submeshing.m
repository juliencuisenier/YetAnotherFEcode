function [Submeshes] = Submeshing(Mesh)

%Give as an output a cell containing Mesh structs of the substructures such
%that : Submeshes{iSub,1} = array containing the coordinates of the nodes
%       Submeshes{iSub,2} = elements array
%       Submeshes{iSub,3} = array containing the original indices of the
%                           global mesh

Submeshes = {};

nSubs = max(Mesh.ELE_TAGS(:,1))-min(Mesh.ELE_TAGS(:,1))+1;

min_index = min(Mesh.ELE_TAGS(:,1)) - 1;

for iSub=1:nSubs
    
    %Find the indices of elements belonging to the i-th substructure
    iElementsIndices = Mesh.ELE_TAGS(:,1)==min_index+iSub;
    
    %Taking the elements
    iElements = Mesh.ELE_NODES(iElementsIndices,:);
    
    %Taking the indices of the nodes of the elements
    iIndices = unique(iElements);
    
    %Reindex the Elements array
    iElements = reindexing_elements_from_global(iElements,iIndices);
    
    %Taking the coordinates of the nodes which belong to the i-th
    %substructure
    iNodes = Mesh.POS(iIndices,:);
    
    %
    iNset = find_nsets(iNodes);
    
    
    %Filling the Submeshes cell
    Submeshes{iSub,1} = iNodes;
    
    Submeshes{iSub,2} = iElements;
    
    Submeshes{iSub,3} = iNset;
    
    Submeshes{iSub,4} = iIndices;
    
end

end
