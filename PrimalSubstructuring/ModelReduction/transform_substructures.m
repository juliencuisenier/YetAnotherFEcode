function transform_substructures(Substructuring)
            % transform_substructures: translates and rotates all
            % substructures to the connection of the master substructure
            % (first colums in interfaces) the master substructure defines
            % the global coordinate system
            tol = 1e-6; % tolerance for the positions
            
            figure('units','normalized','position',[.2 .1 .6 .8])
            hold on %plot of the master substructure
            PlotMesh(Substructuring.Substructures(1).Mesh.nodes,Substructuring.Substructures(1).Mesh.Elements);
            drawnow
            
            for j = 1:(Substructuring.nSubs-1) %go trough all substructures and check for interfaces
            NrSub1 = j;
            for i = NrSub1+1:Substructuring.nSubs
            NrSub2 = i;
            nodepositions1 = Substructuring.Substructures(NrSub1).Mesh.nodes;
            nodepositions2 = Substructuring.Substructures(NrSub2).Mesh.nodes;
            
            
            NodeInt1 = Substructuring.Interfaces(1:end,NrSub1); %extract all possible connections
            NodeInt2 = Substructuring.Interfaces(1:end,NrSub2);
            
            I1 = NodeInt1 ~= 0 ;
            I2 = NodeInt2 ~= 0 ;
            connection = I1.*I2; %boolean array indicating connection between Sub1 and 2
            connection = logical(connection);
            
            NodeInt1 = Substructuring.Interfaces(connection,NrSub1);
            NodeInt2 = Substructuring.Interfaces(connection,NrSub2);
            
            if isempty(NodeInt1)
                continue    %if ther are no connections continue
            end
            originNode1 = NodeInt1(1);  %define the first node as origin for the rotation and translation
            originNode2 = NodeInt2(1);  %(could be random)
            
            origin1 = nodepositions1(originNode1,:);
            origin2 = nodepositions2(originNode2,:);
            
            %traslation vector
            translat = origin1 - origin2;
            %translatation of the second substructure to origin of the
            %first substruct.
            nodepositions2shifted = nodepositions2 + translat;
            Substructuring.Substructures(NrSub2).Mesh.Nodes = nodepositions2shifted;
            
            refNode1 = NodeInt1(end); %take a second node from the interface(could be random)
            refNode2 = NodeInt2(end);
                
            ref1 = nodepositions1(refNode1,:);
            ref2 = nodepositions2shifted(refNode2,:);
            
            if (norm(ref1-ref2))>tol %if they are not yet connected

                vect1 = ref1-origin1;
                vect2 = ref2-origin1;
                
                rot_dir = cross(vect1,vect2)/norm(cross(vect1,vect2));
                rot_angle = atan2(norm(cross(vect1,vect2)),dot(vect1,vect2));
                if abs(abs(vect1/norm(vect1)*(vect2/norm(vect2))')-1)<tol %if the vectors are aligned and the cross product does not exist
                    rot_dir = [1 0 0];
                    rot_angle = pi;
                end
                %compute the angle in this plane
                vect2rot = nodepositions2shifted -origin1; %compute vector from origin to all points
                v_rot = rodrigues_rot(vect2rot,rot_dir,rot_angle); %rodrigues formula function (rotate around normal axis)
                nodepositions2_rotated = v_rot + origin1; %add origin to rotated vectors to have global points
                
                if (norm(ref1-nodepositions2_rotated(refNode2,:)))>tol % atan2 is only defined between
                    v_rot = rodrigues_rot(vect2rot,rot_dir,-rot_angle); % 0-pi and positive
                    nodepositions2_rotated = v_rot + origin1;           % if positions don't match --> rotate by negative angle
                end
               
                Substructuring.Substructures(NrSub2).Mesh.nodes = nodepositions2_rotated;
                
                posInt1 = nodepositions1(NodeInt1,:);
                posInt2 = nodepositions2_rotated(NodeInt2,:);
            else %the first rotation wasn't needed (ref1 and ref2 were already connected)
                posInt1 = nodepositions1(NodeInt1,:);
                posInt2 = nodepositions2shifted(NodeInt2,:);
                nodepositions2_rotated = nodepositions2shifted;
            end
            
            if norm(posInt1-posInt2)>tol % if there is  a position that does not match (normaly true)
                diff = posInt1-posInt2;
                [~,B_ind]=max(diff(:,1).^2+diff(:,2).^2+diff(:,3).^2); % take point with the biggest difference

                Brefnode1 = NodeInt1(B_ind);
                Brefnode2 = NodeInt2(B_ind);
                Bref1 = nodepositions1(Brefnode1,:); 
                Bref2 = nodepositions2_rotated(Brefnode2,:);

                rot_dir = (ref1-origin1)/norm(ref1-origin1); % use the vector between origin and ref as rotation_axis

                Bref1_vect = Bref1- origin1; %generate vector from origin for the B ref points
                Bref2_vect = Bref2- origin1;
                Bref1_vect_ortho = - cross(rot_dir,cross(rot_dir,Bref1_vect)); %only take orthogonal part to the rot-axis
                Bref2_vect_ortho = - cross(rot_dir,cross(rot_dir,Bref2_vect));

                rot_angle = atan2(norm(cross(Bref1_vect_ortho,Bref2_vect_ortho)),...
                    dot(Bref1_vect_ortho,Bref2_vect_ortho)); %compute angle between Bref points (same as before)

                vect2rot = nodepositions2_rotated -origin1;
                v_rot = rodrigues_rot(vect2rot,rot_dir,rot_angle); %rodrigues formula function
                nodepositions2_Brotated = v_rot + origin1; %add origin to rotated vectors

                if (norm(Bref1-nodepositions2_Brotated(Brefnode2,:)))>tol
                v_rot = rodrigues_rot(vect2rot,rot_dir,-rot_angle);
                nodepositions2_Brotated = v_rot + origin1;
                end

                Substructuring.Substructures(NrSub2).Mesh.nodes = nodepositions2_Brotated; %replace node positions with new one

                nodeplot(nodepositions1,Brefnode1,'r' )
                nodeplot(nodepositions2_Brotated,Brefnode2,'r' )
                drawnow
            end
                
            
            nodeplot(nodepositions1,originNode1,'b' )
            nodeplot(nodepositions1,refNode1,'m' )

            PlotMesh(Substructuring.Substructures(NrSub2).Mesh.nodes,Substructuring.Substructures(NrSub2).Mesh.Elements);
            nodeplot(Substructuring.Substructures(NrSub2).Mesh.nodes,originNode2,'b' )
            nodeplot(Substructuring.Substructures(NrSub2).Mesh.nodes,refNode2,'m' )
            drawnow
            end 
            end
        end