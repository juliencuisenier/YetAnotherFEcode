function [Substructures] = find_nsets(Substructures)

%Find the nsets of each substructures in order use them in the Assembly
%class. The nsets will be added in Substructures{iSub,4}

nSubs = size(Substructures,1);

for iSub=1:nSubs
    
    nset_iSub = {};
    
    for j =1:3
        
    nset_iSub{j} = find(Substructures{iSub,1}(:,j)==...
        min(Substructures{iSub,1}(:,j)));
    
    nset_iSub{3+j} = find(Substructures{iSub,1}(:,j)==...
        max(Substructures{iSub,1}(:,j)));
    
    end
    
    Substructures{iSub,4} = nset_iSub;
    
end


end

