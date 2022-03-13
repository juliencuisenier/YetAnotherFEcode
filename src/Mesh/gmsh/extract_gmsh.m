function [nodes,elements,nset] = extract_gmsh(Mesh)

%Extract the informations from the load_gmsh output in order to use yafec
%material

elements = Mesh.ELE_NODES;
    
nodes = Mesh.POS;

nset = {};

for i=1:3
    
    nset{i} = find(nodes(:,i) == min(nodes(:,i)));

    nset{3+i} = find(nodes(:,i) == max(nodes(:,i)));

end
    
end

