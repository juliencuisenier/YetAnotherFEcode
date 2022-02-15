function Substructuring = create_Interface(Substructuring,NrSub1,NrSub2,Inodes1,...
                Inodes2,Bnode1,Bnode2,B2node1,B2node2)% add a connection angle for trishells and beams 
            % Interface: computes  the connecting nodes of two surfaces and
            % stores them in the Interfaces matrix:
            % NrSub1 = number of the first substructure (bigger surface area)
            % NrSub2 = number of the second substructure
            % Inodes1 = nodes of surface 1
            % Inodes2 = nodes of surface 2
            % Bnode 1 = edge or boundary node of the first surface together
            % with Bnode2 of the second surface defines the exact location 
            % of the connection, this two point always match
            %
            % B2node 1 & 2 = additional boundary nodes to better define the
            % interface surface, e.g when symetric
            %
            %In general substructure 2 is the smaller surface and the whole area is 
            %connected. Sub 1 has
            %a bigger connection surface (has not to hold for number of nodes)
            
            nodepositions1 = Substructuring.Substructures(NrSub1).Mesh.nodes;
            nodepositions2 = Substructuring.Substructures(NrSub2).Mesh.nodes;
            
            positionsInt1 = nodepositions1(Inodes1,:);
            positionsInt2 = nodepositions2(Inodes2,:);
            
            startpoint1 = nodepositions1(Bnode1,:);
            startpoint2 = nodepositions2(Bnode2,:);
            % vector from the edged/startpoint to all surface points
            vector1_st = positionsInt1-startpoint1;
            vector1_st = abs(sqrt(vector1_st(:,1).^2+ vector1_st(:,2).^2 + vector1_st(:,3).^2));
            % vector from the edged/startpoint to all surface points
            vector2_st = positionsInt2-startpoint2;
            vector2_st = abs(sqrt(vector2_st(:,1).^2+ vector2_st(:,2).^2 + vector2_st(:,3).^2));
            tol = 1e-6;
            
            if nargin == 7
            
            [vector1_st_sort,vectind1] = sort(vector1_st,'descend');
            [vector2_st_sort,vectind2] = sort(vector2_st,'descend');
            
            foundmeshpoint = 0;
            count_vect2 = 1;
            % extract a unique second node based on the distance radius so
            % that the point matches on both surfaces. With this point a
            % localization of all other points on the surface is possible
            % chose the point as far away from the first edge node to avoid
            % symetries
            while  foundmeshpoint == 0
                count_vect2 = count_vect2 + 1;
                if vector2_st_sort(count_vect2-1)-vector2_st_sort(count_vect2)> tol &&  vector2_st_sort(count_vect2)-vector2_st_sort(count_vect2+1)> tol 
                   
                    for i = 2:(length(vector1_st)-1)
                        if abs(vector2_st_sort(count_vect2)-vector1_st_sort(i))< tol 
                            if vector1_st_sort(i-1)-vector1_st_sort(i)> tol &&  vector1_st_sort(i)-vector1_st_sort(i+1)> tol 
                                foundmeshpoint = 1;
                                count_vect1 = i;
                            end
                        end
                    end  
                end
            end
            
            meshnode1 = Inodes1(vectind1(count_vect1));
            meshnode2 = Inodes2(vectind2(count_vect2));            
            meshpoint1 = positionsInt1(vectind1(count_vect1),:);
            meshpoint2 = positionsInt2(vectind2(count_vect2),:);
            end
            
            if nargin == 9 % used the second localization point from the input
                meshnode1 = B2node1;
                meshnode2 = B2node2;

                meshpoint1 = nodepositions1(B2node1,:);
                meshpoint2 = nodepositions2(B2node2,:);            
            
            end

            % checking all surface nodes if it is a connecting node
            vector1_ms = positionsInt1-meshpoint1;
            vector1_ms = abs(sqrt(vector1_ms(:,1).^2+ vector1_ms(:,2).^2 + vector1_ms(:,3).^2));
            
            vector2_ms = positionsInt2-meshpoint2;
            vector2_ms = abs(sqrt(vector2_ms(:,1).^2+ vector2_ms(:,2).^2 + vector2_ms(:,3).^2));
            
            connectedNodes = 1;
            for i = 1: length(vector1_ms)
                for j = 1:length(vector2_ms)
                    if abs(vector1_st(i)-vector2_st(j))< tol && abs(vector1_ms(i)-vector2_ms(j))< tol
                          local_interface(connectedNodes,1) = Inodes1(i); %store the connected nodes in
                          local_interface(connectedNodes,2) = Inodes2(j); % the local interface matrix
                          connectedNodes = connectedNodes +1;               
                    end
                    
                    
                end
            end
            
%             figure
%             hold on
%             PlotMesh(nodepositions1,Substructuring.Substructures(NrSub1).Mesh.Elements);
%             nodeplot(nodepositions1,Inodes1,'k' )
%             nodeplot(nodepositions1,local_interface(:,1),'m' )
%             nodeplot(nodepositions1,Bnode1,'b')
%             nodeplot(nodepositions1,meshnode1,'g')
%             title(['Output ss#' num2str(NrSub1) ' of "create\_Interface(\_)"'])
%             
%             figure
%             hold on
%             PlotMesh(nodepositions2,Substructuring.Substructures(NrSub2).Mesh.Elements);
%             nodeplot(nodepositions2,Inodes2,'k' )
%             nodeplot(nodepositions2,local_interface(:,2),'m' )
%             nodeplot(nodepositions2,Bnode2,'b')
%             nodeplot(nodepositions2,meshnode2,'g')
%             title(['Output ss#' num2str(NrSub2) ' of "create\_Interface(\_)"'])
            
            %store the new interface nodes in the general Interfaces matrix
            newInterface = zeros(size(local_interface,1),Substructuring.nSubs);
            newInterface(:,NrSub1) = local_interface(:,1);
            newInterface(:,NrSub2) = local_interface(:,2);
            Substructuring.Interfaces = [Substructuring.Interfaces; newInterface];
            Substructuring.nInt = Substructuring.nInt();
        end