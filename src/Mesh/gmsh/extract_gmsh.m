function [nodes,elements,nset] = extract_gmsh(Mesh)

%Extract the informations from the load_gmsh output in order to use yafec
%material

elements = Mesh.HEXAS(:,1:8);

indices = unique(elements);
    
nodes = Mesh.POS(indices,:);

nset = {};

for i=1:3
    
nset{i} = find(nodes(:,i) == min(nodes(:,i)));

nset{3+i} = find(nodes(:,i) == max(nodes(:,i)));

end


%Reindexing elements -> gmsh take account of 1D and 2D elements in the
%count of indices. We have to reindexed the element array in order to make
%it consistent with the indices of the nodes array

[nbElm, nbNodes_by_Elm] = size(elements);

reindexed_elements = zeros(nbElm, nbNodes_by_Elm);

for jElm=1:nbElm
    
    for kNode=1:nbNodes_by_Elm
        new_indice = find(indices == elements(jElm,kNode));
        
        reindexed_elements(jElm,kNode)=new_indice;     
    end
        
end

elements = reindexed_elements;
    
end

