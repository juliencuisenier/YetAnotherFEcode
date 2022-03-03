function [Substructures] = find_nsets(Substructures)

%Find the nsets of each substructures in order use them in the Assembly
%class. The nsets will be added in Substructures{iSub,4}

nSubs = size(Substructures,1);

for iSub=1:nSubs
    
    nset_iSub = {};
    
    nset_iSub{1} = find(Substructures{iSub,1}(:,1)==...
        min(Substructures{iSub,1}(:,1)));
    
    nset_iSub{2} = find(Substructures{iSub,1}(:,2)==...
        min(Substructures{iSub,1}(:,2)));
    
    nset_iSub{3} = find(Substructures{iSub,1}(:,3)==...
        min(Substructures{iSub,1}(:,3)));
    
    nset_iSub{4} = find(Substructures{iSub,1}(:,1)==...
        max(Substructures{iSub,1}(:,1)));
    
    nset_iSub{5} = find(Substructures{iSub,1}(:,2)==...
        max(Substructures{iSub,1}(:,2)));
    
    nset_iSub{6} = find(Substructures{iSub,1}(:,3)==...
        max(Substructures{iSub,1}(:,3)));
    
    Substructures{iSub,4} = nset_iSub;
    
end


end

