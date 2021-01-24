% curvefit_simulation_toDSC_trigfunction.m
%
% Nathanael Kazmierczak, 12/25/2019

NUMDISKS = 12;
AVERAGE_FLAG = 0;
FIT_FLAG = 1;

if ~exist('stack4D_output','var')
    addpath('/Volumes/NPKResearch/BediakoLab/NathanaelLargeMatlab');
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/Code_from_Colin/BlinkingSimulations/bilayergraphenesimulations/NPK');
    filename = 'reconTBLG1p2Moire_r-ang0p8total_r-d40.mat';
    load(filename);
    filename2 = '12252019diskcenters_for_rastersimulation_reconTBLG1p2Moire_r-ang0p8total_r-d40.mat';
    load(filename2);
    filename3 = '12252019_simulation_disk_averages.mat';
    load(filename3);
end

if ~exist('disk_centers','var')
    [ disk_centers, radius ] = makeDiskCenters( stack4D_output(:,:,1,1) );
end

plotDP(stack4D_output(:,:,1,1));
viscircles(disk_centers,radius*ones(size(disk_centers,1),1));

[xs,ys,r,c] = size(stack4D_output);

if AVERAGE_FLAG
    % Actually pull down the disk averages data that we wanted.
    xbase = 1:xs;
    ybase = 1:ys;
    disk_average_storage = zeros(r,c,NUMDISKS);
    [xspace,yspace] = meshgrid(xbase,ybase);
    for rr = 1:r
        fprintf('Getting disk averages of data row %d out of %d...\n',rr,r);
        for cc = 1:c
            thisDP = double(stack4D_output(:,:,rr,cc));
            %         global_row_index = (i-1)*chunksize + rr;
            for q = 1:size(disk_centers,1)
                x0 = disk_centers(q,1);
                y0 = disk_centers(q,2);
                tf = isInCircle(xspace,yspace,x0,y0,radius);
                disk_average_storage(rr,cc,q) = mean(thisDP(tf));
            end
        end
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
    DSC_fit_storage = zeros(r,c,2);
    residuals_storage = zeros(r,c,12);
    RMSR_storage = zeros(r,c);
    % normalize and zero
    mins = min(min(disk_average_storage));
    disk_average_storage = disk_average_storage - mins;
    maxes = repmat(max(max(disk_average_storage)),[r,c,1]);
    disk_average_storage = disk_average_storage./maxes;
    %     weight_vector = [1,1,0,1,1,0,...
    %         1,1,1,1,1,1];   % Because I don't particularly like the outlier-esque look of inner disks 3 and 6.
    weight_vector = [1,1,1,1,1,1,...
        1,1,1,1,1,1];
    %     weight_vector = [1,1,1,1,1,1,...
    %         0,0,0,0,0,0];
    %     weight_vector = [0,0,0,0,0,0,...
    %         1,1,1,1,1,1];
    
    for i = 1:r
        fprintf('Optimizing data row %d out of %d...\n',i,r);
        for j = 1:c
            this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
            % Crudely, start with the [0,0] guess. May need to multistart
            % this if convergence seems suspect.
            initialDSC_guess = [0,0];
            % lsqnonlin syntax -- computing the vector values of the
            % objective function and not the RMSR.
            fit_scaling = 0.9;
            objfcn = @(DSCvector) (this_disk_averages - trigFittingFunctions(DSCvector,fit_scaling)).*weight_vector;
%             objfcn = @(DSCvector) (this_disk_averages - blinkingFittingFunction(DSCvector)).*weight_vector;

%             lb = [-1.5,-1.5];
            lb = [-1.5,0];   % This bound means no negative y-values, to break the optimization ambiguity.
            ub = [1.5,1.5];
            %             DSC_fit_result = lsqnonlin(objfcn,initialDSC_guess,lb,ub,options);
            [ DSC_fit_result ] = multistartDiskOptimization( objfcn, lb, ub, options );
            DSC_fit_storage(i,j,:) = permute(DSC_fit_result,[3,1,2]);
            this_residuals = objfcn(DSC_fit_result);
            residuals_storage(i,j,:) = permute(this_residuals,[3,1,2]);
            RMSR_storage(i,j) = rms(objfcn(DSC_fit_result));
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
    DSCangle = atan2(DSC2,DSC1);
    
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
%     figure; imagesc(DSCangle);
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
    
    
end
