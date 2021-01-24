classdef DisplacementEquivalenceClass < handle
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        reduced_zone_displacement
        equivalent_displacements
        returnval
    end
    
    methods
        function obj = DisplacementEquivalenceClass(init_displacement,returnval,permutation_number,xbounds,ybounds)
            obj.reduced_zone_displacement = init_displacement;
%             perm_val = 2;
            equivalent_points = DisplacementEquivalenceClass.getEquivalentDisplacements(init_displacement,permutation_number);
            c = equivalent_points(:,1) < xbounds(2) | equivalent_points(:,1) > xbounds(1) | equivalent_points(:,2) < ybounds(2) | equivalent_points(:,2) > ybounds(1);
            equivalent_points(c,:) = [];
            obj.equivalent_displacements = equivalent_points;
            obj.returnval = returnval;
        end
        
        function difference = minus(obj,otherobj)
            [r,c] = size(otherobj.equivalent_displacements);
            dm = repmat(obj.reduced_zone_displacement,r,1);
            differences = dm - otherobj.equivalent_displacements;
            distances = sum((differences).^2,2).^0.5;
            [~,idx] = min(distances);
            candidate_diff = differences(idx,:);
            % The part where we need some convention to break the ambiguity
            % on which way the vector is pointing.
            if candidate_diff(1) < 0
                candidate_diff = -candidate_diff;
            end
            switch obj.returnval
                case 'x'
                    difference = candidate_diff(1);
                case 'y'
                    difference = candidate_diff(2);
                otherwise
                    difference = candidate_diff;
            end
        end
            
        function vecnorm = reportReducedZoneMagnitude(obj)
            vecnorm = norm(obj.reduced_zone_displacement);
        end
        
        
        function showDisplacement(obj)
            if isempty(obj.equivalent_displacements)
                perm_val = 1;
                obj.getEquivalentDisplacements(perm_val);
            end
            figure
            axh = gca;
            plot(obj.equivalent_displacements(:,1),obj.equivalent_displacements(:,2),'o','MarkerEdgeColor',[1,0,0]);
            hold on;
            plot(obj.reduced_zone_displacement(1),obj.reduced_zone_displacement(2),'o','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
            plotFullDisplacementHexagons(axh);
            xlim([-4,4]);
            ylim([-4,4]);
        end
        
        
    end
    
    
    
    
    methods(Static)
        function equivalent_points = getEquivalentDisplacements(reduced_zone_displacement,perm_val)
            basis = permn(-perm_val:perm_val,2);
            % In this function, I think we do want to include the reduced
            % zone displacement in the all-points representation.
%             basis(basis(:,1) == 0 & basis(:,2) == 0,:) = [];
            [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
            w1 = v1 + v2;
            w2 = 2*v2 - v1;
            W = [w1',w2'];
            genpoints = W*basis';
            genpoints_plus = genpoints + repmat(reduced_zone_displacement',1,size(genpoints,2));
            genpoints_minus = genpoints - repmat(reduced_zone_displacement',1,size(genpoints,2));
            equivalent_points = [genpoints_plus,genpoints_minus,-reduced_zone_displacement']';
        end
    end
    
end

