classdef DisplacementEquivalenceClassMagnet < DisplacementEquivalenceClass
    %UNTITLED12 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        indices
        current_displacement
    end
    
    methods
        function obj = DisplacementEquivalenceClassMagnet(init_displacement,indices,permutation_number,xbounds,ybounds)
            obj = obj@DisplacementEquivalenceClass(init_displacement,'all',permutation_number,xbounds,ybounds);
            obj.indices = indices;
            obj.current_displacement = init_displacement;
        end
        
        function H = evaluateCurrentEnergy(obj,displacement_field,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling)
            H = obj.evaluateEnergyAtVector(obj.current_displacement,displacement_field,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
        end
        
        function [best_energy,best_vector] = transitionToBestNewVector(obj,displacement_field,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling,transition_flag)
            numvecs = size(obj.equivalent_displacements,1);
            energies = zeros(numvecs,1);
            for j = 1:numvecs
                this_candidate_vector = obj.equivalent_displacements(j,:);
                H = obj.evaluateEnergyAtVector(this_candidate_vector,displacement_field,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
                energies(j) = H;
            end
            [best_energy,idx] = min(energies);
            best_vector = obj.equivalent_displacements(idx,:);
            if transition_flag
                obj.current_displacement = best_vector;
            end
        end
        
        function H = evaluateEnergyAtVector(obj,vector_of_interest,displacement_field,n_coupling,fixed_emitters_displacements,fixed_emitters_indices,f_coupling)
            [r,c,h] = size(displacement_field);
            H = 0;
            if n_coupling > 0
            if obj.indices(1) == 1 || obj.indices(1) == r || obj.indices(2) == 1 || obj.indices(2) == c
                % Then there will be no neighbor energy interactions.
            else
                neighbors = displacement_field(obj.indices(1)-1:obj.indices(1)+1,obj.indices(2)-1:obj.indices(2)+1,:);
                neighbors_x = neighbors(:,:,1);
                neighbors_y = neighbors(:,:,2);
                neighbors_linear = [neighbors_x(:), neighbors_y(:)];
%                 neighbor_mags = sqrt(neighbors_x.^2 + neighbors_y.^2);
%                 this_mag = sqrt(vector_of_interest(1).^2 + vector_of_interest(2).^2);
%                 n_for_dot = [neighbors_x(:),neighbors_y(:)];
                % weight by vector size
%                 n_maxes = max(neighbor_mags,repmat(this_mag,9,1),[],2);
%                 n_mins = min(neighbor_mags,repmat(this_mag,9,1),[],2);
%                 scalar_weights = n_mins ./ n_maxes;  % The point is to penalize when the vectors are different from each other.
%                 dot_vals = n_for_dot*vector_of_interest';  % gives the directional component.
                this_ds_distance = sum((neighbors_linear-repmat(vector_of_interest,9,1)).^2,2).^0.5;
                H = H + sum(this_ds_distance*n_coupling);
            end
            end
            % Now evaluate energy from the fixed point emitters
            
            if f_coupling > 0
%             num_emitters = size(fixed_emitters_displacements,1);
%             disp('begin an emitter calc');
%             for i = 1:num_emitters
                this_rs_distance = sum((fixed_emitters_indices-obj.indices).^2,2).^0.5;
                this_displacement_distance = sum((fixed_emitters_displacements-vector_of_interest).^2,2).^0.5;
%                 directionality = dot(obj.indices,fixed_emitters_directions(i,:));
                emitter_energies = f_coupling(1)*this_rs_distance.^(-f_coupling(2)).*this_displacement_distance;
                H = H + sum(emitter_energies);
%             end
            end
%             H
        end
            
        
    end
    
end

