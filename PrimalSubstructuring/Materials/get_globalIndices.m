function [globalIndices] = get_globalIndices(Submeshes)
%GET_GLOBALINDICES 
%gives as an output a cell containing the global indices of each submesh.

nSubs = size(Submeshes,1);

globalIndices = cell(1,nSubs);

for iSub=1:nSubs
    
    globalIndices{iSub} = Submeshes{iSub,4};

end

