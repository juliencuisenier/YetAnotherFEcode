function Substructuring = create_Interface(Substructuring,Submeshes,Sub1, Sub2)

    %Define the interface between Sub1 and Sub2
    %OUTPUTS : Substructuring.Interface and Substructuring.nInt
    %INPUTS : Substructuring -> Substructuring struct 
    %         Submeshes -> cell containing the the information of the 
    %                      substructure's mesh 
    %         Sub1 -> index of one of the substucture
    %         Sub2 -> index of the other substructure
    
    indices1 = Submeshes{Sub1,3};
    
    indices2 = Submeshes{Sub2,3};
    
    Interface = intersect(indices1,indices2);
    
    reindexed_Interface1 = [];
    reindexed_Interface2 = [];
    
    for i=1:size(Interface)
       reindexed_Interface1 = [reindexed_Interface1; ...
           find(Submeshes{Sub1,3}== Interface(i))];
    end
    
    for i=1:size(Interface)
       reindexed_Interface2 = [reindexed_Interface2; ...
           find(Submeshes{Sub2,3}== Interface(i))];
    end
    
    reindexed_Interface = [reindexed_Interface1 reindexed_Interface2];
    
    Substructuring.Interfaces = [Substructuring.Interfaces reindexed_Interface];
    
    Substructuring.nInt = Substructuring.nInt();

end