function indices_corrector = corr_gmsh_indices(PrimalSub,GlobalMesh)
% corr_gmsh_indices makes a vector that make correspondance the index of the
% substructuring method and the index of the gmsh method.
% Let i an index of a global vector of the substructuring class, and make
% it localize in a vector us of a substructure. Now i is in the j-th
% position of  us. We can know thanks to globalIndices where the j-th index
% is in the global mesh, thus i.
%
%In the corr_gmsh_indices vector we have :
% [...]
% [index of a substructuring global vector] <-in the index which correspond
% [...]                                      to the right index in the
%                                            global gmsh mesh



us = PrimalSub.Us;
indices_corrector = zeros(PrimalSub.nDOFglobal,1);

for iSub=1:PrimalSub.nSubs
    for j=1:PrimalSub.Substructures(iSub).Mesh.nNodes       
        gmsh_index = PrimalSub.globalIndices{iSub}(j);
        gmsh_dofs = GlobalMesh.get_DOF_from_nodeIDs(gmsh_index);
        indices_corrector(gmsh_dofs) = us{iSub}(3*(j-1)+1:3*(j-1)+3);
    end
end

end