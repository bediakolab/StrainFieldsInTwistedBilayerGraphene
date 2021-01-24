function [ final_reddispfield ] = shiftFittingBoundsMagnet( init_reddispfield, use_init, neighbor_range, nomod_maxint, distance_flag )
% shiftFittingBounds

%% Prepare for shifting logic
if nomod_flag
    if nargin > 4
        MAXINT = nomod_maxint;
    else
        MAXINT = 10;
    end
else
    MAXINT = 4; % For the reduced shifting, this works. For nomod reconstruction, will undoubtedly need more.
end

if isempty(donotuseroi)
    donotuseroi = zeros(size(init_reddispfield,1),size(init_reddispfield,2));
end

basis = permn(-MAXINT:MAXINT,2);
[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
w1 = v1 + v2;
w2 = 2*v2 - v1;
W = [w1',w2'];
genpoints = W*basis';

%% Build the fitting unit cell
vec1 = vectorLine([1,0],[0,0]);
vec2 = vectorLine([1,0],v2+v1);
vec3 = vectorLine(2*v1-v2,[0,0]);
vec4 = vectorLine(2*v1-v2,v1+v2);

%%
if ~isempty(AAregions)
    AAdists = sum((AAregions - AAregions(1,:)).^2).^(1/2);
    [AAminval] = min(AAdists);
    AAregions = sortrows(AAregions);
    disp('hi');
end

%% Loop over vectors
[r,c,h] = size(init_reddispfield);
if use_init
    computed = true(r,c);
    final_reddispfield = init_reddispfield;
else
    computed = false(r,c);
    final_reddispfield = zeros(r,c,h);
end
for i = 1:r
    fprintf('Processing displacement field row %d of %d.\n',i,r);
    for j = 1:c
        if donotuseroi(i,j)
            continue
        end
        thisvec = permute(init_reddispfield(i,j,:),[3,1,2]);
        % Evaluate the options in the fitting unit cell to see which
        % satisfy the geometry
        genpoints_plus = genpoints + repmat(thisvec,1,size(genpoints,2));
        genpoints_minus = genpoints - repmat(thisvec,1,size(genpoints,2));
        genpoints_all = [genpoints_plus,genpoints_minus,-thisvec]';
        
        if ~nomod_flag
            c1 = vec1 < genpoints_all;
            c2 = vec2 > genpoints_all;
            c3 = vec3 < genpoints_all;
            c4 = vec4 > genpoints_all;
            %         idxs = find();
            candidate_new_vecs = genpoints_all(c1 & c2 & c3 & c4,:);
            num_candidates = size(candidate_new_vecs,1);
        elseif distance_flag == 5
            mc1 = i > AAregions(:,1);
            AAregions_red = AAregions;
            AAregions_red(~mc1,:) = [];
            mc2 = j > AAregions(:,2);
            idx = find(mc2,1);
            anchorpt = AAregions_red(idx,:)
            
            vec1 = vectorLine([1,0],[0,0]+anchorpt);
            vec2 = vectorLine([1,0],v2+v1+anchorpt);
            vec3 = vectorLine(2*v1-v2,[0,0]+anchorpt);
            vec4 = vectorLine(2*v1-v2,v1+v2+anchorpt);
            
            c1 = vec1 < genpoints_all;
            c2 = vec2 > genpoints_all;
            c3 = vec3 < genpoints_all;
            c4 = vec4 > genpoints_all;
            %         idxs = find();
            candidate_new_vecs = genpoints_all(c1 & c2 & c3 & c4,:);
            num_candidates = size(candidate_new_vecs,1);
        else
            candidate_new_vecs = genpoints_all;
            num_candidates = size(candidate_new_vecs,1);
        end
        
        figure
        plot(candidate_new_vecs(:,1),candidate_new_vecs(:,2),'ro');
        axh = gca;
        plotFullDisplacementHexagons(axh);
        
        % Key step: iterate through the candidate vectors to ascertain
        % which will be the best fit for the given location in the raster.
        % First, obtain all of the neighboring vectors.
        idxvec = -neighbor_range:neighbor_range;
        pixelidxs = permn(idxvec,2);
        pixelidxs(pixelidxs(:,1) == 0 & pixelidxs(:,2) == 0,:) = [];
        theseinds = [i,j];
        newinds = repmat(theseinds,size(pixelidxs,1),1) + pixelidxs;
        newinds(newinds(:,1) <= 0 | newinds(:,2) <= 0 | newinds(:,1) > r | newinds(:,2) > c,:) = [];
        % First, check if the neighbors have been computed yet.
        lininds = sub2ind([r,c],newinds(:,1),newinds(:,2));
        newinds(~computed(lininds),:) = [];
%         newinds2 = newinds;
%         newinds = [newinds,ones(size(newinds,1),1)];
%         newinds2 = [newinds2,2*ones(size(newinds2,1),1)];
%         newinds = [newinds;newinds2];
%         lininds = sub2ind([r,c,h],newinds(:,1),newinds(:,2),newinds(:,3));
        neighbors = zeros(0,2);
        for z = 1:size(newinds,1)
            neighbors(end+1,1) = final_reddispfield(newinds(z,1),newinds(z,2),1);
            neighbors(end,2) = final_reddispfield(newinds(z,1),newinds(z,2),2);
        end
        num_neighbors = size(neighbors,1);
        if isempty(neighbors)
            num_neighbors = 0;
        end
        
        distances = zeros(num_candidates,1);
        lastinds = theseinds - [1,1];
        if any(lastinds < 1)
            lastvec = [0,0];
        else
            lv1 = final_reddispfield(1:lastinds(1),1:lastinds(2),1);
            lv2 = final_reddispfield(1:lastinds(1),1:lastinds(2),2);
            lastvec = [lv1(:),lv2(:)];
        end
        
        for k = 1:num_candidates
            this_candidate = candidate_new_vecs(k,:);
            if distance_flag == 4
                if sum(this_candidate.^2)^(1/2) <= mean(sum(lastvec.^2).^(1/2))
                    distances(k) = 1000;
                    continue
                end
            else
                distances(k) = 0;
            end
            
            
            if distance_flag == 6
                if ~use_init
                    i1 = max(1,i-neighbor_range);
                    i2 = max(1,j-neighbor_range);
                    %                 i1 = 1;
                    %                 i2 = 1;
                    test_grad_mat_x = final_reddispfield(i1:i,i2:j,1);
                    test_grad_mat_y = final_reddispfield(i1:i,i2:j,2);
                    insertion = size(test_grad_mat_x);
                    test_grad_mat_x(insertion(1),insertion(2)) = this_candidate(1);
                    test_grad_mat_y(insertion(1),insertion(2)) = this_candidate(2);
                else
%                     istartbase = i-neighbor_range;
%                     iendbase = i+neighbor_range;
%                     jstartbase = j-neighbor_range;
%                     jendbase = j+neighbor_range;
                    
                    final_reddispfield_temp = final_reddispfield;
                    final_reddispfield_temp(i,j,1) = this_candidate(1);
                    final_reddispfield_temp(i,j,2) = this_candidate(2);
                    
                    istart = max(1,i-neighbor_range);
                    iend = min(r,i+neighbor_range);
                    jstart = max(1,j-neighbor_range);
                    jend = min(c,j+neighbor_range);
                    %                 i1 = 1;
                    %                 i2 = 1;
                    test_grad_mat_x = final_reddispfield_temp(istart:iend,jstart:jend,1);
                    test_grad_mat_y = final_reddispfield_temp(istart:iend,jstart:jend,2);
%                     insertion = [istart,mean([jstart,jend])];
%                     test_grad_mat_x(insertion(1),insertion(2)) = this_candidate(1);
%                     test_grad_mat_y(insertion(1),insertion(2)) = this_candidate(2);
                end
                    
                    
                if size(test_grad_mat_x,1) > 1 && size(test_grad_mat_x,2) > 1
                    [xx,xy] = gradient(test_grad_mat_x);
                    [yx,yy] = gradient(test_grad_mat_y);
                    gradnorm = mean([rms(rms(xx))^2,rms(rms(xy))^2,rms(rms(yx))^2,rms(rms(yy))^2]);
                else
                    xderiv = gradient(test_grad_mat_x);
                    yderiv = gradient(test_grad_mat_y);
                    gradnorm = mean([rms(rms(xderiv)),rms(rms(yderiv))]);
                end
                distances(k) = gradnorm;
                
                    
            else
                for q = 1:num_neighbors
                    this_neighbor = neighbors(q,:);
                    if this_neighbor(1) == 0 && this_neighbor(2) == 0
                        continue
                    end
                    if distance_flag == 1 || distance_flag == 5
                        thisdist = sum((this_neighbor - this_candidate).^2).^(1/2);
                        %                     if i == 17 && j == 19
                        %                         disp('hi');
                        %                     end
                    elseif distance_flag == 2 || distance_flag == 4
                        this_candidate_amp = sum(this_candidate.^2).^(1/2);
                        this_neighbor_amp = sum(this_neighbor.^2).^(1/2);
                        thisdist = sum((this_neighbor_amp - this_candidate_amp).^2).^(1/2);
                        
                    elseif distance_flag == 3
                        this_candidate_angle = atan2(this_candidate(2),this_candidate(1));
                        this_neighbor_angle = atan2(this_neighbor(2),this_neighbor(1));
                        thisdist = sin(abs(this_candidate_angle - this_neighbor_angle));
                    elseif distance_flag == 7
                        this_candidate_amp = sum(this_candidate.^2).^(1/2);
                        this_neighbor_amp = sum(this_neighbor.^2).^(1/2);
                        thisdistamp = sum((this_neighbor_amp - this_candidate_amp).^2).^(1/2);
                        
                        this_candidate_angle = atan2(this_candidate(2),this_candidate(1));
                        this_neighbor_angle = atan2(this_neighbor(2),this_neighbor(1));
                        thisdistangle = sin(abs(this_candidate_angle - this_neighbor_angle));
                        
                        thisdist = 0.5*thisdistamp + 0.5*thisdistangle;
                    end
                    
                    
                    distances(k) = distances(k) + thisdist;
                end
            end
        end
        
        if i == 42 && j == 21
            disp('hi');
        end
        
        if i == 1 && j == 1
            [~,bestcandidateidx] = min(sum(candidate_new_vecs.^2,2).^(1/2));
        else
            [~,bestcandidateidx] = min(distances);
        end
        chosenvec = candidate_new_vecs(bestcandidateidx,:);
        final_reddispfield(i,j,:) = permute(chosenvec,[2,3,1]);
        computed(i,j) = true;
    end
end

end

