function [reindexed_elements] = reindexing_elements_from_global(Elements,Indices)

%Change the indices of the nodes in the elements array, in order to 
%correspond them to the new indexing caused by the submeshing

%OUTPUTS : reindex_Elements -> Array containing the new relevant indices
%INPUTS : Elements -> Elements' array
%         Indices -> vector containing the global indices


[nbElm, nbNodes_by_Elm] = size(Elements);

reindexed_elements = zeros(nbElm, nbNodes_by_Elm);

for iElm=1:nbElm
    
       for jNode=1:nbNodes_by_Elm
        
          new_indice = find(Indices == Elements(iElm,jNode));  
          reindexed_elements(iElm,jNode) = new_indice;    
            
       end    
       
end
    
    
    
        
    
end




