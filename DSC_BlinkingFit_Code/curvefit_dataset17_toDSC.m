% curvefit_dataset17_toDSC.m
%
% Nathanael Kazmierczak, Dec 2019

AVERAGE_DISKS_FLAG = 0;
FIT_FLAG = 1;

NUMDISKS = 12;
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/New2019_1023/Dataset1');
filename = '1_MV_1p5_9_40x40_ss0p5nm_1s_spot9_alpha=1_bin2_cl=130_60kV_normalized.h5';

count = [1024,1024,4,40];  % stays the same while start varies.
fullsize = [1024,1024,40,40];
NUMCHUNKS = 10;
chunksize = fullsize(3)/NUMCHUNKS;

if AVERAGE_DISKS_FLAG
    load('12232019Dataset1DiskCenters.mat');
    
    disk_average_storage = zeros(fullsize(3),fullsize(4),NUMDISKS);
    for i = 1:NUMCHUNKS
        fprintf('Beginning chunk %d.\n',i);
        thisstart = [1,1,(i-1)*chunksize+1,1];
        thisdata = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',thisstart,count);
        
        
        if i == 1 && ~exist('disk_centers','var')
            while true
                
                [ disk_centers, radius ] = makeDiskCenters( thisdata );
                randDPs = randi(3,3,2);
                plotDP(thisdata(:,:,randDPs(1,1),randDPs(1,2)));
                viscircles(disk_centers,radius*ones(size(disk_centers,1),1));
                plotDP(thisdata(:,:,randDPs(2,1),randDPs(2,2)));
                viscircles(disk_centers,radius*ones(size(disk_centers,1),1));
                plotDP(thisdata(:,:,randDPs(3,1),randDPs(3,2)));
                viscircles(disk_centers,radius*ones(size(disk_centers,1),1));
                
                tf = input('1/0: Are the defined disks satisfactory? ');
                if tf
                    input('Save these disk definitions.');
                    break
                else
                    disp('Clearing disk center choices.');
                    
                    clear disk_centers radius
                end
                
            end
        end
        
        % Actually pull down the disk averages data that we wanted.
        [~,~,r,c] = size(thisdata);
        xbase = 1:fullsize(1);
        ybase = 1:fullsize(2);
        [xspace,yspace] = meshgrid(xbase,ybase);
        for rr = 1:r
            for cc = 1:c
                thisDP = thisdata(:,:,rr,cc);
                global_row_index = (i-1)*chunksize + rr;
                for q = 1:size(disk_centers,1)
                    x0 = disk_centers(q,1);
                    y0 = disk_centers(q,2);
                    tf = isInCircle(xspace,yspace,x0,y0,radius);
                    disk_average_storage(global_row_index,cc,q) = mean(thisDP(tf));
                end
            end
        end
        clear thisdata
    end
    
    makeDataBlinkingPlot( disk_average_storage(:,:,1) )
    title 'Inner disk 1'
    makeDataBlinkingPlot( disk_average_storage(:,:,2) )
    title 'Inner disk 2'
    makeDataBlinkingPlot( disk_average_storage(:,:,3) )
    title 'Inner disk 3'
    makeDataBlinkingPlot( disk_average_storage(:,:,4) )
    title 'Inner disk 4'
    makeDataBlinkingPlot( disk_average_storage(:,:,5) )
    title 'Inner disk 5'
    makeDataBlinkingPlot( disk_average_storage(:,:,6) )
    title 'Inner disk 6'
    makeDataBlinkingPlot( disk_average_storage(:,:,7) )
    title 'Outer disk 1'
    makeDataBlinkingPlot( disk_average_storage(:,:,8) )
    title 'Outer disk 2'
    makeDataBlinkingPlot( disk_average_storage(:,:,9) )
    title 'Outer disk 3'
    makeDataBlinkingPlot( disk_average_storage(:,:,10) )
    title 'Outer disk 4'
    makeDataBlinkingPlot( disk_average_storage(:,:,11) )
    title 'Outer disk 5'
    makeDataBlinkingPlot( disk_average_storage(:,:,12) )
    title 'Outer disk 6'
    
end


if FIT_FLAG
    % Perform the fit. Can also be done separately.
    options = optimoptions('lsqnonlin');
    options.Display = 'off';
    diskavname = '12232019Dataset1DiskAverages_Normalized.mat';
    load(diskavname);
    DSC_fit_storage = zeros(40,40,2);
    residuals_storage = zeros(40,40,12);
    RMSR_storage = zeros(40,40);
    % normalize and zero
    mins = min(min(disk_average_storage));
    disk_average_storage = disk_average_storage - mins; 
    maxes = repmat(max(max(disk_average_storage)),[40,40,1]);
    disk_average_storage = disk_average_storage./maxes;
%     weight_vector = [1,1,0,1,1,0,...
%         1,1,1,1,1,1];   % Because I don't particularly like the outlier-esque look of inner disks 3 and 6.
    weight_vector = [1,1,1,1,1,1,...
        1,1,1,1,1,1];
%     weight_vector = [1,1,1,1,1,1,...
%         0,0,0,0,0,0];
%     weight_vector = [0,0,0,0,0,0,...
%         1,1,1,1,1,1];

    for i = 1:fullsize(3)
        fprintf('Optimizing data row %d out of %d...\n',i,fullsize(3));
        for j = 1:fullsize(4)
            this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
            % Crudely, start with the [0,0] guess. May need to multistart
            % this if convergence seems suspect.
            initialDSC_guess = [0,0];
            % lsqnonlin syntax -- computing the vector values of the
            % objective function and not the RMSR.
            objfcn = @(DSCvector) (this_disk_averages - blinkingFittingFunction(DSCvector)).*weight_vector;
            
            lb = [-1.5,-1.5];
            ub = [1.5,1.5];
%             DSC_fit_result = lsqnonlin(objfcn,initialDSC_guess,lb,ub,options);
            [ DSC_fit_result ] = multistartDiskOptimization( objfcn, lb, ub, options );
            DSC_fit_storage(i,j,:) = permute(DSC_fit_result,[3,1,2]);
            this_residuals = objfcn(DSC_fit_result);
            residuals_storage(i,j,:) = permute(this_residuals,[3,1,2]);
            RMSR_storage(i,j) = rms(objfcn(DSC_fit_result));
        end
    end
end

% Visualize the range of values converged to.
figure 
DSC1 = DSC_fit_storage(:,:,1);
DSC2 = DSC_fit_storage(:,:,2);
scatter(DSC1(:),DSC2(:),'Filled');
xlabel('Cartesian DSC component 1');
ylabel('Cartesian DSC component 2');
title('Converged DSC scatterplot');

% Make spatial DSC plots (Code taken from BilayerGraphene -- maybe should
% refactor at some point?)
DSCamp = (DSC1.^2 + DSC2.^2).^0.5;
DSCangle = atan2(DSC1,DSC2);

figh = figure; contourf(DSCamp,50,'LineStyle','None');
colormap(fire)
axis equal
shading interp
colorbar
title(sprintf(strcat(['DSC field amplitude'])));
% set(gcf, 'Position',  [0, 0, 500, 800])
xlabel('x');
ylabel('y');

figure; contourf(DSCangle,50,'LineStyle','None');
colormap(hsv)
axis equal
shading interp
colorbar
title(sprintf(strcat(['DSC field angle'])));
% set(gcf, 'Position',  [0, 0, 500, 800])
xlabel('x');
ylabel('y');

% Visualize residuals
figure; 
surfc(RMSR_storage);
colormap(jet);
colorbar
xlabel('x');
ylabel('y');
title('Fit RMSR');
