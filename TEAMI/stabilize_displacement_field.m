% stabilize_displacement_field.m
%
% Attempt to remove displacement field rotational problems for plotting and
% strain mapping. Using dataset 4 from the TEAM I session, day 1.
%
% Nathanael Kazmierczak, 02/29/2020

% load('TEAMI_Day1_DS7_objectdata.mat');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset7RadialResidualsInterp');
load('objectdataRadialResiduals.mat');
dfield = m4.DSC_fit_storage;
subdivnum = 80;
dfield = dfield(1:subdivnum,1:subdivnum,1:2);

PLOT_INIT = 1;
if PLOT_INIT
    damp = (dfield(:,:,1).^2 + dfield(:,:,2).^2).^(1/2);
    dangle = atan2(dfield(:,:,2),dfield(:,:,1));
    figh1 = figure;
    imagesc(damp); colormap(fire); axis equal; colorbar;
    title('Initial amp');
    figure;
    imagesc(dangle); colormap(hsv); axis equal; colorbar;
    title('Initial angle');
end

% use_init = 1;
% nomod_maxint = 5;
% neighbor_range = 5;
% nomod_flag = 1;
% distance_flag = 3
TEST_PICKS = false;
if TEST_PICKS
    

    AAregions = zeros(0,2);
    while true
        disp('Please click an AA region.');
        AAregions(end+1,:) = ginput(1);
        choice = input('Are you done clicking AA regions? 1/0');
        if choice
            break
        end
    end
      
else
    AAregions = [];
end


ROI = false;
if ROI
    figure(figh1);
    disp('Select the ignoring region.');
    donotuseroi = roipoly;
end

BOUNDARY_MAP = false;
if BOUNDARY_MAP
    Hsize = 7;
    ampthresh = 1;
    medfilt_damp = medfilt2(damp);
    sub_amp = abs(damp - medfilt_damp) > ampthresh;
    % Do not permit changes around the edges
    sub_amp(1,:) = 0;
    sub_amp(end,:) = 0;
    sub_amp(:,1) = 0;
    sub_amp(:,end) = 0;
    damp(sub_amp) = medfilt_damp(sub_amp);
    figure;
    imagesc(damp); colormap(fire); axis equal; colorbar;
    title('NPK Median filter amplitudes');
    
    H = ones(Hsize);
    H = H/sum(sum(H));
%     smoothed_damp = filter2(H,damp);
    smoothed_damp = imgaussfilt(damp,8);
%     smoothed_damp(1:ceil(Hsize)/2,:) = damp(1:ceil(Hsize)/2,:);
    figh2 = figure;
    [ damp_trimmed, rowindices, colindices ] = trimArray( smoothed_damp,4 );
    imagesc(damp_trimmed); colormap(fire); axis equal; colorbar;
    title('3x3 moving average amplitudes after medfilt2');
    set(gca,'ydir','normal');
    
    
    
%     H = ones(3);
%     smoothed_damp = imgaussfilt(H,medfilt_damp);
%     figh2 = figure;
%     imagesc(smoothed_damp); colormap(fire); axis equal; colorbar;
%     title('3x3 moving average amplitudes after medfilt2');

    while true
        disp('Click a location on the graph to start a minimum path-following point.');
        [x,y] = ginput(1);
        BW_init = false(size(smoothed_damp));
        [ BW_final ] = followMinPath( damp_trimmed, BW_init, [x,y] );
        hold on
        spy(BW_final);
    end
end


MAGNET = true;
if MAGNET
    % Build the cell array of DisplacementMagnets;
    [r,c,h] = size(dfield);
    rbase = 1:r;
    cbase = 1:c;
    [rgrid,cgrid] = meshgrid(rbase,cbase);
    magnets = cell(r,c);
    energies = zeros(r,c);
    for i = 1:subdivnum
        for j = 1:subdivnum
            magnets{i,j} = DisplacementEquivalenceClassMagnet(permute(dfield(i,j,:),[1,3,2]),[i,j]);
        end
    end
    % Set the fixed emitter points.close
    fighs = m4.makeDisplacementMapsFromBlinking(0,[],0,'flat');
    close(fighs(3:end));
    figure(fighs(2));
    [c,d] = ginput(6);
    fixed_emitters_indices = [d,c];  % This needs to be swtiched because we will use it to compare distances against indices
    disp('Please click on the displacement vectors for each of the six emitters, in order.');
    figure(fighs(1));
    [a,b] = ginput(6);
    fixed_emitters_directions = [a,b];
    
    %%
    
    figure
    dfieldx = dfield(:,:,1);
    dfieldy = dfield(:,:,2);
    quiver(rgrid(:),cgrid(:),dfieldx(:),dfieldy(:),0);
    title(sprintf('Starting quiver'));
    
    for i = 1:subdivnum
        for j = 1:subdivnum
            magnets{i,j} = DisplacementEquivalenceClassMagnet(permute(dfield(i,j,:),[1,3,2]),[i,j]);
        end
    end
    
    % Compute starting energies
    n_coupling = 1;
    f_coupling = [100000,0.5];
    for i = 1:subdivnum
        for j = 1:subdivnum
            energies(i,j) = magnets{i,j}.evaluateCurrentEnergy(dfield,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
        end
    end
    figure
    imagesc(energies); colormap(fire);
    starting_energy = sum(sum(energies));
    colorbar;
    title(sprintf('Beginning Energy, energy %f',starting_energy));
    
    niter = 5;
    transition_flag = 1;
%     figure
    for q = 1:niter
        q
        if q > 1
            f_coupling = [10000,1];
            n_coupling = 10;
        end
        if q > 2
            f_coupling = [10,1];
            n_coupling = 100;
        end
            
        for i = 1:subdivnum
            for j = 1:subdivnum
%                 energies(i,j) = magnets{i,j}.evaluateCurrentEnergy(n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
                [energies(i,j),~] = magnets{i,j}.transitionToBestNewVector(dfield,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling,transition_flag);
            end
        end
        figure
        imagesc(energies); colormap(fire);
        this_energy = sum(sum(energies));
        hold on;
        colorbar;
        title(sprintf('Iteration %d, energy %f',q,this_energy));
        
        current_dfield = zeros(80,80,2);
        for i = 1:80
            for j = 1:80
                current_dfield(i,j,:) = permute(magnets{i,j}.current_displacement,[1,3,2]);
            end
        end
        figure
        dfieldx = current_dfield(:,:,1);
        dfieldy = current_dfield(:,:,2);
        quiver(rgrid(:),cgrid(:),dfieldx(:),dfieldy(:),0);
        hold on 
        quiver(fixed_emitters_indices(:,2),fixed_emitters_indices(:,1),fixed_emitters_directions(:,1),fixed_emitters_directions(:,2),0,'r');
        title(sprintf('Iteration %d',q));
        dfield = current_dfield;
        
    end
end





% nomod_flag = 1;
% use_init = 0;
% nomod_maxint = 3;
% neighbor_range = 3;
% distance_flag = 1;
% [ processed_dfield ] = shiftFittingBounds( dfield, use_init, neighbor_range, nomod_flag, nomod_maxint, distance_flag, AAregions, donotuseroi );
% [ processed_dfield ] = shiftFittingBounds( flipud(dfield), use_init, neighbor_range, nomod_flag );
% [ processed_dfield ] = shiftFittingBounds( processed_dfield, use_init, neighbor_range, nomod_flag, nomod_maxint, distance_flag  );
% processed_dfield = flipud(processed_dfield);

% use_init = 1;
% distance_flag = 2;
% neighbor_range = 5;
% [ processed_dfield ] = shiftFittingBounds( processed_dfield, use_init, neighbor_range, nomod_flag, nomod_maxint, distance_flag, AAregions );




