function [Substructures] = reindexing_elements_from_global(Substructures)

%Change the indices of the nodes in the elements array, in order to 
%correspond them to the new indexing caused by the submeshing

%OUTPUTS : Substructures -> Cell with the {iSub,2} updated with local 
%                           indices
%INPUTS : Substructures -> Cell containing the informations of each 
%                          substructure

nSubs = size(Substructures,1);

for iSub=1:nSubs
    
    iElements = Substructures{iSub,2};
    
    iIndices = Substructures{iSub,3};

    
    [nbElm, nbNodes_by_Elm] = size(iElements);

    reindexed_elements = zeros(nbElm, nbNodes_by_Elm);

    for jElm=1:nbElm
    
        for kNode=1:nbNodes_by_Elm
        
            if iElements(jElm,kNode)==0
                break;
            end
        
        new_indice = find(iIndices == iElements(jElm,kNode));
        
        reindexed_elements(jElm,kNode)=new_indice;
            
        end
        
    end
    
    Substructures{iSub,2} = reindexed_elements;
    
        
end
    
end




