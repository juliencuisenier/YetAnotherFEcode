function [nset] = find_nsets(Nodes)

%Find the nsets of the nodes array in order use them in the Assembly
%class. 

nset = cell(1,6);
    
for j =1:3
        
    nset{j} = find(Nodes(:,j)==...
        min(Nodes(:,j)));
    
    nset{3+j} = find(Nodes(:,j)==...
        max(Nodes(:,j)));
    
end


end

