classdef ReducedAssembly < Assembly
    properties
        % Already contains properties in the Assembly class
        
        V       % Reduction basis
    end

    methods
        function self = ReducedAssembly(Mesh,V)
            % Assembly (superclass) constructor
            self@Assembly(Mesh);
            
            % set reduction basis
            self.V = V;
        end

        function set.V(self,V)
            if size(V,1) ~= self.Mesh.nDOFs %#ok<*MCSUP>
                warning(['Reduction basis size incorrect: should have ' ...
                    num2str(self.Mesh.nDOFs) ' rows'])
            end
            self.V = V;
        end

        function [f] = uniform_body_force(self,varargin)
            m = size(self.V,2);
            f = zeros(m,1);
            V = self.V;             %#ok<*PROP>
            [elementWeights,~] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % extracting domain elements from the supplied set
            domainElements = ~[self.Mesh.Elements(elementSet).isBoundary];
            
            Elements = self.Mesh.Elements;
            parfor j = elementSet(domainElements)
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                fe = thisElement.uniformBodyForce;
                f = f + elementWeights(j) * (Ve.' * fe);
            end 
        end
        
     

        function [K] = matrix(self,elementMethodName,varargin)
            % This function assembles a generic finite element matrix from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            K = zeros(m,m);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                Ke = thisElement.(elementMethodName)(inputs{:});
                K = K + elementWeights(j) * (Ve.' * Ke * Ve);
            end
        end

        function [f] = vector(self,elementMethodName,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system

            m = size(self.V,2);
            f = zeros(m,1);

            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            parfor j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                fe = thisElement.(elementMethodName)(inputs{:});
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end
        
        function [K, f] = matrix_and_vector(self,elementMethodName, varargin)
            % This function assembles a generic finite element matrix and
            % vector from its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            
            K = zeros(m,m);
            f = zeros(m,1);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions
            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                [Ke, fe] = thisElement.(elementMethodName)(inputs{:});
                K = K + elementWeights(j) * (Ve.' * Ke * Ve);
                f = f + elementWeights(j) * (Ve.' * fe);
            end
        end
    
        function [T] = tensor(self,elementMethodName,SIZE,sumDIMS,varargin)
            % This function assembles a generic finite element vector from
            % its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level vector Fe.
            % For this to work, a method named elementMethodName which
            % returns the appropriate vector must be defined for all
            % element types in the FE Mesh.            
            % NOTE: in this function, we reduce all the dimensions with 
            % the same reduction basis
            
            m = size(self.V,2);
            [~,I] = find(SIZE == m);
            T = tenzeros(SIZE);
            Elements = self.Mesh.Elements;
            V = self.V;                %#ok<*PROPLC>
            
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            % Computing element level contributions

            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:); %#ok<*PFBNS>
                [Te, ~] = thisElement.(elementMethodName)(inputs{:});                
                % transform tensor
                 Vcell = cell(ndims(Te),1);
                 Vcell(:) = {Ve.'};
                 T = T + elementWeights(j) * ttm(Te, Vcell,I);
%                  Vt = tensor(Ve');    % reduction basis (transpose)
%                             Vnt  = tensor(Ve);     % reduction basis
%             if elementMethodName=='tensor_Q3'
%                 T  = T + elementWeights(j)*ttt(ttt(ttt(Vt,Te ,2,1),Vnt,3,1),Vnt ,2,1);
%             elseif elementMethodName=='tensor_Q4'
%                 T= T + elementWeights(j)*ttt(ttt(ttt(ttt(Vt,Te  ,2,1),Vnt,4,1),Vnt ,3,1),Vnt ,2,1);
%             end
            
            [subs, T] = sparsify(T,[],sumDIMS);
            T = sptensor(subs, T, SIZE);
        end

        end
        
            function [G, b] = constructGbass(self,qq, varargin)
            % This function assembles a generic finite element matrix and
            % vector from its element level counterpart.
            % elementMethodName is a string input containing the name of
            % the method that returns the element level matrix Ke.
            % For this to work, a method named elementMethodName which
            % returns the appropriate matrix must be defined for all
            % element types in the FE Mesh.            
            % NOTE: it is assumed that the input arguments are provided in
            % the full (unreduced) system
            
            m = size(self.V,2);
            
            K = zeros(m,m);
            f = zeros(m,1);
            Elements = self.Mesh.Elements;
            V = self.V;
            % parsing element weights
            [elementWeights,inputs] = self.parse_inputs(varargin{:});
            
            % extracting elements with nonzero weights
            elementSet = find(elementWeights);
            
            nt=size(qq,2);
            G=sparse(zeros(m*nt,self.Mesh.nElements));
            b=sparse(zeros(m*nt,1));
            % Computing element level contributions
            for k=1:nt
             
                z2=k*m;
                z1=z2-m;
            for j = elementSet
                thisElement = Elements(j).Object;
                index = thisElement.iDOFs;          
                Ve = V(index,:);
                [Ke, fe] = thisElement.('tangent_stiffness_and_force')(V*qq(:,k));
                
               
                f = elementWeights(j) * (Ve.' * fe);
                G(z1+1:z2,elementSet(j))=f;
                 G=sparse(G);
            end
           
           b(z1+1:z2)=sum( G(z1+1:z2,:),2);
           b=sparse(b);
            end
        1
            end
    end
end
