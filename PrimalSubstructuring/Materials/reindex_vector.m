function [reindexed_v] = reindex_vector(corr_gmsh_indices,v)
%REINDEX_VECTOR interchanges the values of v in order to make them accurate
%with the global gmsh mesh indexation. Function corr_gmsh_indices must has
%been computed before using this function.

nDOFs = length(corr_gmsh_indices);

reindexed_v = zeros(nDOFs,1);

for i=1:nDOFs
    reindexed_v(i)=v(corr_gmsh_indices(i));
end

end

