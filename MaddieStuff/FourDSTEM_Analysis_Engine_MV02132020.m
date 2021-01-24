classdef FourDSTEM_Analysis_Engine_MV02132020 < handle
    % Class for performing the blinking analysis of a 4DSTEM dataset
    %
    % Nathanael Kazmierczak, 01/01/2019
    
    properties
        filename
        defined_disk_filename
        saved_plots_folderpath
        saved_plots_foldername
        probe_filename
        datacube_size
        num_load_chunks
        CHUNK_LOAD_WIDTH
        disk_centers
        disk_radius
        disk_averages
        skipped_disks
        residuals_storage
        RMSR_storage
        DSC_fit_storage
        scan_stepsize  % in units of nm, which is what is given by Maddie's filenames.
        exx  % strain map components
        eyy
        exy
        eyx
        gxy
        theta_torsion
        xaxis
        yaxis
        power_law_fit_mask
        hBN_mask
        beamstop_mask
        graphene_mask
        detector_mask
        power_law_fit_parameters
        optimal_power_law_fit
        disk_centers_G1
        disk_centers_G2
        indices_G1
        indices_G2
        lattice_vectors_G1
        lattice_vectors_G2
        beam_center_coords
        sFit_struct
        weight_vector
        average_BN_probe
    end
    
    methods
        % in the constructor, resultsfolderpath gives the full path to the
        % folder where images will be saved, including the name of the
        % folder itself.
        function obj = FourDSTEM_Analysis_Engine_MV02132020(filename,disk_filename,resultsfolderpath,scan_stepsize,probe_filename)
            % Assume that the filepath has been set in the driver
            % disk_filename is optional
            obj.filename = filename;
            if nargin > 1
                obj.defined_disk_filename = disk_filename;
            end
            if nargin > 3
                obj.scan_stepsize = scan_stepsize;  % The input is in nanometers, keep it that way
            end
            if nargin > 4
                obj.probe_filename = probe_filename;
            end
            
            try
                h = h5info(filename);
                h1 = h.Groups.Groups.Groups;  % Because of the weird bit with multiple outputs.
                h2 = h1.Groups;
                h3 = h2.Datasets.Dataspace;
                obj.datacube_size = h3.Size;
            catch
                warning('File is not a valid .h5. Partial loading will be disabled.');
            end
            
            % determine number of load chunks
            obj.CHUNK_LOAD_WIDTH = 5;
            obj.num_load_chunks = ceil(obj.datacube_size(3)/obj.CHUNK_LOAD_WIDTH);
            
            splitted = strsplit(resultsfolderpath,'/');
            numsplit = length(splitted);
            thispath = '';
            for i = 1:(numsplit-1)
                thispath = strcat([thispath,splitted{i},'/']);
            end
            obj.saved_plots_folderpath = thispath;
            %             obj.saved_plots_folderpath = resultsfolderpath;
            obj.saved_plots_foldername = splitted{end};
            
            % Set x and y axes
            obj.xaxis = obj.scan_stepsize*(0:(obj.datacube_size(4)-1));
            obj.yaxis = obj.scan_stepsize*(0:(obj.datacube_size(3)-1));
        end
        
        
        
        % Define the integration regions for the blinking calculations.
        function [] = setBlinkingDisks(obj,explicit_chunks)
            disp('Beginning the setBlinkingDisks() method...');
            if nargin < 2
                explicit_chunks = [];
            end
            
            if ~isempty(obj.defined_disk_filename)
                fprintf('It appears disks for this file have already been saved at %s.\n',obj.defined_disk_filename);
                tf = input('Do you wish to (0) load or (1) overwrite these disks?');
                if ~tf
                    disp('Loading disks');
                    load(obj.defined_disk_filename);
                    obj.disk_centers = disk_centers;
                    obj.disk_radius = radius;
                    obj.skipped_disks = skipped_disks;
                    return
                end
            end
            
            if ~isempty(explicit_chunks)
                first_load = explicit_chunks(1);
                middle_load = explicit_chunks(2);
            else
                first_load = 1;
                middle_load = round(obj.num_load_chunks/2);
            end
            [data1,~,~] = obj.partialLoad(first_load);
            [data2,~,~] = obj.partialLoad(middle_load);
            diffraction_patterns = cat(3,data1,data2);
            %             [ obj.disk_centers, obj.disk_radius ] = makeDiskCenters( diffraction_patterns );
            
            while true
                [ disk_centers_temp, radius, skipped_disks ] = makeDiskCenters( diffraction_patterns );
                randDPs_dim3 = randi(obj.CHUNK_LOAD_WIDTH,3,1);
                randDPs_dim4 = randi(obj.datacube_size(4),3,1);
                plotDP(diffraction_patterns(:,:,randDPs_dim3(1),randDPs_dim4(1)));
                viscircles(disk_centers_temp,radius*ones(size(disk_centers_temp,1),1));
                plotDP(diffraction_patterns(:,:,randDPs_dim3(2),randDPs_dim4(2)));
                viscircles(disk_centers_temp,radius*ones(size(disk_centers_temp,1),1));
                plotDP(diffraction_patterns(:,:,randDPs_dim3(3),randDPs_dim4(3)));
                viscircles(disk_centers_temp,radius*ones(size(disk_centers_temp,1),1));
                
                tf = input('1/0: Are the defined disks satisfactory? ');
                if tf
                    disp('Saving these disk definitions under the following filename:');
                    s = strsplit(obj.filename,'.');
                    obj.defined_disk_filename = sprintf('%s__integration_disks.mat',s{1});
                    fprintf('%s\n',obj.defined_disk_filename);
                    disk_centers = disk_centers_temp;
                    save(obj.defined_disk_filename,'disk_centers','radius','skipped_disks');
                    break
                else
                    disp('Clearing disk center choices.');
                    clear disk_centers_temp radius
                end
            end
            obj.disk_centers = disk_centers_temp;
            obj.disk_radius = radius;
            obj.skipped_disks = skipped_disks;
        end
        
        
        
        function plotIntegrationDisks(obj)
            [diffraction_patterns,~,~] = obj.partialLoad(randi(obj.num_load_chunks));
            randDPs_dim3 = randi(obj.CHUNK_LOAD_WIDTH,3,1);
            randDPs_dim4 = randi(obj.datacube_size(4),3,1);
            plotDP(diffraction_patterns(:,:,randDPs_dim3(1),randDPs_dim4(1)));
            viscircles(obj.disk_centers,obj.disk_radius*ones(size(obj.disk_centers,1),1));
            plotDP(diffraction_patterns(:,:,randDPs_dim3(2),randDPs_dim4(2)));
            viscircles(obj.disk_centers,obj.disk_radius*ones(size(obj.disk_centers,1),1));
            plotDP(diffraction_patterns(:,:,randDPs_dim3(3),randDPs_dim4(3)));
            viscircles(obj.disk_centers,obj.disk_radius*ones(size(obj.disk_centers,1),1));
        end
        
        
        
        function disk_average_storage = integrateDisks(obj,subtract_power_law_flag)
            disp('Beginning the integrateDisks() method...');
            NUMDISKS = 12;
            disk_average_storage = zeros(obj.datacube_size(3),obj.datacube_size(4),NUMDISKS);
            for i = 1:obj.num_load_chunks
                [thisdata,indices_start,indices_end] = obj.partialLoad(i);
                % Actually pull down the disk averages data that we wanted.
                [~,~,r,c] = size(thisdata);
                xbase = 1:obj.datacube_size(1);
                ybase = 1:obj.datacube_size(2);
                [xspace,yspace] = meshgrid(xbase,ybase);
                for rr = 1:r
                    for cc = 1:c
                        thisDP = thisdata(:,:,rr,cc);
                        if subtract_power_law_flag
                            thisDP = double(thisDP) - obj.optimal_power_law_fit;
                        end
                        global_row_index = indices_start(3) + rr - 1;
                        for q = 1:size(obj.disk_centers,1)
                            x0 = obj.disk_centers(q,1);
                            y0 = obj.disk_centers(q,2);
                            tf = isInCircle(xspace,yspace,x0,y0,obj.disk_radius);
                            disk_average_storage(global_row_index,cc,q) = mean(thisDP(tf));
                        end
                    end
                end
                obj.disk_averages = disk_average_storage;
                clear thisdata
            end
        end
        
        
        
        function makeBlinkingPlots(obj,disk_ids,saveplotflag,shadingtype)
            if nargin < 2 || isempty(disk_ids)
                figh = makeMultiDataBlinkingPlot( obj.disk_averages(:,:,[1:3,7:9]),shadingtype,obj.xaxis,obj.yaxis);
            else
                figh = makeMultiDataBlinkingPlot( obj.disk_averages(:,:,disk_ids),shadingtype,obj.xaxis,obj.yaxis);
            end
            set(gcf, 'Position',  [0, 0, 1400, 800])
            if saveplotflag
                currentd = pwd;
                %cd(obj.saved_plots_folderpath);
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                saveas(figh,'BlinkingPlots.png');
                cd(currentd);
            end
        end
        
        
        
        % Throw in options here about the weighting vector, the fitting
        % function, the scaling factor, etc.
        function fitBlinking(obj,trig_prefactor,weight_vector_override)
            
            disp('Beginning the fitBlinking() method...');
            if nargin < 2
                trig_prefactor = 0.9;
            end
            if nargin < 3
                weight_vector_override = [];
            end
            
            options = optimoptions('lsqnonlin');
            options.Display = 'off';
            %             diskavname = '12232019Dataset1DiskAverages_Normalized.mat';
            %             load(diskavname);
            if isempty(obj.disk_averages)
                error('Please define the disks using obj.setBlinkingDisks()');
            else
                disk_average_storage = obj.disk_averages;
            end
            DSC_fit_storage = zeros(obj.datacube_size(3),obj.datacube_size(4),2);
            residuals_storage = zeros(obj.datacube_size(3),obj.datacube_size(4),12);
            RMSR_storage = zeros(obj.datacube_size(3),obj.datacube_size(4));
            % normalize and zero
            mins = min(min(disk_average_storage));
            % NPK 02/02/2020: test getting rid of the zeroing, as this may
            % be biasing the fit.
            %             disk_average_storage = disk_average_storage - mins;
            maxes = repmat(max(max(disk_average_storage)),[obj.datacube_size(3),obj.datacube_size(4),1]);
            disk_average_storage = disk_average_storage./maxes;
            %     weight_vector = [1,1,0,1,1,0,...
            %         1,1,1,1,1,1];   % Because I don't particularly like the outlier-esque look of inner disks 3 and 6.
            %             weight_vector = [1,1,1,1,1,1,...
            %                 1,1,1,1,1,1];
            if isempty(weight_vector_override)
                weight_vector = ones(1,12);
            else
                weight_vector = weight_vector_override;
            end
            % Not removing the skipped disks will crash the program.
            weight_vector(obj.skipped_disks) = 0;
            weight_vector
            obj.weight_vector = weight_vector;
            
            %     weight_vector = [1,1,1,1,1,1,...
            %         0,0,0,0,0,0];
            %     weight_vector = [0,0,0,0,0,0,...
            %         1,1,1,1,1,1];
            
            %             makeMultiDataBlinkingPlot( disk_average_storage(:,:,[1:3,7:9]));
            
            for i = 1:obj.datacube_size(3)
                fprintf('Optimizing data row %d out of %d...\n',i,obj.datacube_size(3));
                for j = 1:obj.datacube_size(4)
                    this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
                    this_disk_averages(isnan(this_disk_averages)) = 0;  % trust the weighting vector has taken care of it.
                    % Crudely, start with the [0,0] guess. May need to multistart
                    % this if convergence seems suspect.
                    %                     initialDSC_guess = [0,0];
                    % lsqnonlin syntax -- computing the vector values of the
                    % objective function and not the RMSR.
                    %             objfcn = @(DSCvector) (this_disk_averages - blinkingFittingFunction(DSCvector)).*weight_vector;
                    %                     trig_prefactor = 0.9;
                    objfcn = @(DSCvector) (this_disk_averages - trigFittingFunctions(DSCvector,trig_prefactor,0)).*weight_vector;
                    lb = [-1.24,0];
                    trigfun_flag = 1;
                    %             lb = [-1.5,-1.5];
                    
                    ub = [1.24,1.43];  % really a/sqrt(3) and a/2, I believe.
                    %             DSC_fit_result = lsqnonlin(objfcn,initialDSC_guess,lb,ub,options);
                    
                    [ DSC_fit_result ] = multistartDiskOptimization( objfcn, lb, ub, options, trigfun_flag );
                    DSC_fit_storage(i,j,:) = permute(DSC_fit_result,[3,1,2]);
                    this_residuals = objfcn(DSC_fit_result);
                    residuals_storage(i,j,:) = permute(this_residuals,[3,1,2]);
                    RMSR_storage(i,j) = rms(objfcn(DSC_fit_result));
                end
            end
            
            obj.residuals_storage = residuals_storage;
            obj.RMSR_storage = RMSR_storage;
            obj.DSC_fit_storage = DSC_fit_storage;
        end
        
        
        
        % These used to be in with fitBlinking() itself.
        function makeDisplacementMapsFromBlinking(obj,saveplotflag,median_filter_threshold,individual_residual_flag,shading_type)
            % Visualize the range of values converged to.
            if nargin < 3
                median_filter_threshold = 0;
            end
            if nargin < 4
                individual_residual_flag = 0;
            end
            figh_DSCscatter = figure;
            DSC1 = obj.DSC_fit_storage(:,:,1);
            DSC2 = obj.DSC_fit_storage(:,:,2);
            if ~isempty(median_filter_threshold)
                [ DSC1, DSC2 ] = npkMedFilt2Displacement( DSC1, DSC2, median_filter_threshold(1), median_filter_threshold(2) );
            end
            scatter(DSC1(:),DSC2(:),'Filled');
            xlabel('Cartesian DSC component 1');
            ylabel('Cartesian DSC component 2');
            title('Converged DSC scatterplot');
            
            % Make spatial DSC plots (Code taken from BilayerGraphene -- maybe should
            % refactor at some point?)
            DSCamp = (DSC1.^2 + DSC2.^2).^0.5;
            DSCangle = atan2(DSC2,DSC1);
            
            figh_DSCamp = figure;
            if strcmp(shading_type,'interp') 
                contourf(DSCamp,50,'LineStyle','None');
            else
                imagesc(DSCamp);
            end
            colormap(fire)
            axis equal
            shading interp
            colorbar
            title(sprintf(strcat(['DSC field amplitude (angstroms)'])));
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x');
            ylabel('y');
            set(gca,'ydir','normal');
            
            %             figure; contourf(DSCangle,50,'LineStyle','None');
            figh_DSCangle = figure; 
            if strcmp(shading_type,'interp')
                contourf(DSCangle,50,'LineStyle','None');
            else
                imagesc(DSCangle);
            end
            colormap(hsv)
            axis equal
            shading interp
            colorbar
            title(sprintf(strcat(['DSC field angle (radians)'])));
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x');
            ylabel('y');
            set(gca,'ydir','normal');
            
            % Visualize residuals
            RMSR_storage = obj.RMSR_storage;
            figh_DSCresiduals = figure;
            surfc(RMSR_storage);
            colormap(jet);
            colorbar
            xlabel('x');
            ylabel('y');
            title('Fit RMSR');
            
            if individual_residual_flag
                diskres_fighs = zeros(1,12);
                disk_strings = cell(1,12);
                for i = 1:12
                    if obj.weight_vector(i)
                        figh = figure;
                        this_disk_residual = obj.residuals_storage(:,:,i);
                        surfc(this_disk_residual);
                        title(sprintf('DSC Residuals for Disk %d',i));
                        colorbar;
                        diskres_fighs(i) = figh;
                        disk_strings{i} = sprintf('Disk%dResidualsPlot.png',i);
                    end
                end
            end
            
            if saveplotflag
                currentd = pwd;
                %cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                saveas(figh_DSCscatter,'DSCscatter.png');
                saveas(figh_DSCamp,'DSCamplitude.png');
                saveas(figh_DSCangle,'DSCangle.png');
                saveas(figh_DSCresiduals,'DSCresiduals.png');
                if individual_residual_flag
                    for i = 1:12
                        if obj.weight_vector(i)
                            saveas(diskres_fighs(i),disk_strings{i});
                        end
                    end
                end
                cd(currentd);
            end
            
        end
        
        
        
        % Obtaining the numerical gradient.
        function makeStrainMapsFromBlinking(obj,saveplotflag,shadingtype,colormaptype,trim_edge_pixel_num,median_filter_threshold)
            if nargin < 3
                median_filter_threshold = 0;
            end
            if isempty(obj.DSC_fit_storage)
                error('Run fitBlinking() before makeStrainMaps().');
            end
            NUM_CONTOUR_LINES = 50;
            % DSC_fit_storage is a 3-way array in which the first component
            % of the 3rd dimension is ux(x,y) and the second is uy(x,y)
            uX = obj.DSC_fit_storage(:,:,1);
            uY = obj.DSC_fit_storage(:,:,2);
            
            if ~isempty(median_filter_threshold)
                [ uX, uY ] = npkMedFilt2Displacement( uX, uY, median_filter_threshold(1), median_filter_threshold(2) );
            end
            
            spacing = obj.scan_stepsize;  % This is in nm and we want to keep it that way.
            uX = 0.1*uX;  % converts the X displacement to nm
            uY = 0.1*uY;  % converts the X displacement to nm
            [obj.exx,obj.exy] = gradient(uX,spacing);
            [obj.eyx,obj.eyy] = gradient(uY,spacing);
            obj.exx = obj.exx - mean(mean(obj.exx));
            obj.eyy = obj.eyy - mean(mean(obj.eyy));
            obj.exy = obj.exy - mean(mean(obj.exy));
            obj.eyx = obj.eyx - mean(mean(obj.eyx));
            obj.gxy = obj.exy + obj.eyx;
            obj.theta_torsion = 0.5*(obj.eyx - obj.exy);
            
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            xbase((end-trim_edge_pixel_num+1):end) = [];
            xbase(1:trim_edge_pixel_num) = [];
            ybase((end-trim_edge_pixel_num+1):end) = [];
            ybase(1:trim_edge_pixel_num) = [];
            
            % Using the following function to trim edges (which may have
            % artifacts)
            %             [ A_trimmed ] = trimArray( A,trim_width )
            figh_exx = figure;
            plot_exx = trimArray(obj.exx,trim_edge_pixel_num);
            
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_exx);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_exx,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            colormap(colormaptype)
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('exx strain component (nm x displacment / nm x real space)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            figh_eyy = figure;
            plot_eyy = trimArray(obj.eyy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_eyy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_eyy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            colormap(colormaptype)
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('eyy strain component (nm y displacment / nm y real space)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            figh_gxy = figure;
            plot_gxy = trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_gxy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_gxy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            colormap(colormaptype)
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('gxy total shear strain (nm x(y) displacment / nm y(x) real space)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            figh_theta = figure;
            plot_theta = trimArray(obj.theta_torsion,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_theta);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_theta,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            colormap(colormaptype)
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('Torsional strain: theta (nm x(y) displacment / nm y(x) real space)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            if saveplotflag
                currentd = pwd;
                %cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                saveas(figh_exx,'exx_strain.png');
                saveas(figh_eyy,'eyy_strain.png');
                saveas(figh_gxy,'gxy_total_shear_strain.png');
                saveas(figh_theta,'theta_torsional_strain.png');
                cd(currentd);
            end
        end
        
        
        
        function makeQuiverPlotFromBlinking(obj,saveplotflag)
            figh = figure;
            DSCx = obj.DSC_fit_storage(:,:,1);
            DSCy = obj.DSC_fit_storage(:,:,2);
            [xspace,yspace] = meshgrid(obj.xaxis,obj.yaxis);
            quiver(xspace(:),yspace(:),DSCx(:),DSCy(:),0,'k');
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Displacement field');
            %             set(gca,'ydir','reverse')
            axis equal
            
            if saveplotflag
                currentd = pwd;
                %cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                saveas(figh,'displacement_quiverplot.png');
                cd(currentd);
            end
        end
        
        
        
        
        function makeInteractiveDisplacementDarkFieldPlot(obj)
            DSC1 = obj.DSC_fit_storage(:,:,1);
            DSC2 = obj.DSC_fit_storage(:,:,2);
            f = figure;
            
            set(gcf, 'Position',  [0, 0, 800, 400])
            s1 = subplot(1,2,1);
            %             ax1 = axes('Parent',s1);
            h1 = scatter(DSC1(:),DSC2(:),'Filled');
            title('Interactive Displacement Plot');
            % axis tight
            ylim([-0.1,1.5]);
            xlim([-1.5,1.5]);
            xlabel('X displacement');
            ylabel('Y displacement');
            
            s2 = subplot(1,2,2);
            %             ax2 = axes('Parent',s2);
            boolmat = false(size(DSC1));
            colormap(gray);
            h2 = imagesc(obj.xaxis,obj.yaxis,boolmat);
            axis equal
            xlim([obj.xaxis(1)-0.5*obj.scan_stepsize, obj.xaxis(end)+0.5*obj.scan_stepsize]);
            ylim([obj.yaxis(1)-0.5*obj.scan_stepsize, obj.yaxis(end)+0.5*obj.scan_stepsize]);
            xlabel('x (nm)');
            ylabel('y (nm)');
            %             set(gca,'position',[0.02 0.04 0.96 0.92]);
            set(s1,'position',[0.1,0.1,0.4,0.8]);
            
            e = imrect(s1);
            fcn = @(pos) obj.darkFieldCallback(pos,s2);
            fcn(e.getPosition());
            id = addNewPositionCallback(e,fcn);
        end
        
        
        % Saves the current instance of the FourDSTEM_Analysis_Engine to
        % the results folder so that you do not have to re-run the script
        % to manipulate the fitting results.
        function saveObject(obj)
            name = inputname(1);
            command = sprintf('save(''objectdata.mat'',''%s'')',name);
            
            currentd = pwd;
            %cd(obj.saved_plots_folderpath)
            if ~exist(obj.saved_plots_foldername,'dir')
                mkdir(obj.saved_plots_foldername);
            end
            cd(obj.saved_plots_foldername);
            evalin('base',command);
            cd(currentd);
        end
        
        
        
        function fitElasticPowerLaw(obj,power_law_mask)
            %             [data,~,~] = obj.partialLoad(1);
            %             data = data(:,:,2,2);  % Sometimes the first one is unusually bright for some reason.
            if nargin < 2
                power_law_mask = [];
            end           
            
            [data] = double(obj.singleLoad(2,2));  % Sometimes the first one is unusually bright for some reason.
            % For the hBN
            if isempty(obj.power_law_fit_mask) && isempty(power_law_mask)
                
                disp('Prepare to mask off the hBN lattice and the beamstop.');
                [ ~, hBN_mask, ~, beamstop_mask, ~ ] = makeHexagonalLatticeMask( data );
                disp('Prepare to mask off the graphene lattice.');
                [ ~, graphene_mask, ~, ~, ~ ] = makeHexagonalLatticeMask( {data} );
                disp('Prepare to demarcate the detector field of view');
                disp('Click once on the detector center and once for the detector radius.');
                [xc,yc] = ginput(2);
                imagesize = size(data);
                rowbase = 1:imagesize(1);
                colbase = 1:imagesize(2);
                [rowspace,colspace] = meshgrid(rowbase,colbase);
                r = sqrt((xc(1)-xc(2))^2 + (yc(1)-yc(2))^2);
                detector_mask = ~isInCircle(colspace,rowspace,yc(1),xc(1),r);
                
                obj.hBN_mask = hBN_mask;
                obj.beamstop_mask = beamstop_mask;
                obj.graphene_mask = graphene_mask;
                obj.detector_mask = detector_mask;
                obj.power_law_fit_mask = hBN_mask | beamstop_mask | graphene_mask | detector_mask;
                
                obj.plotPowerLawMask();
            end
            if ~isempty(power_law_mask)
                input('Preparing to use user-supplied power law mask. Press any key to continue.');
                obj.power_law_fit_mask = power_law_mask;
            end
            
            data(obj.power_law_fit_mask) = -1;
            %             use_teeny_mask = 0;
            %             teeny_mask_radius_factor = 0;
            %             draw_hard_delete_region = 2;
            %             THRESHOLD = 150;
            %             [beamstop_masks,Irng,Jrng,maskeddata] = makeBeamstopCenterMasks(data,use_teeny_mask,teeny_mask_radius_factor,draw_hard_delete_region,THRESHOLD);
            %
            %             [ masked_image, power_law_fit_mask, out_lattice_points, beamstop_mask, lattice_vectors ] = makeHexagonalLatticeMask( maskeddata, radius, init_lattice_points, init_beamstop_mask );
            %
            %             maskeddata = data(
            powerLawFit = @(c) getFullRadialPowerLawFit(c,data);
            y0_init = round(mean([1,size(data,1)]));
            x0_init = round(mean([1,size(data,2)]));
            log10A_init = 15;
            B_init = -4;
            c_initial_guesses = [x0_init,y0_init,log10A_init,B_init];
            options = optimset;
            options.Display = 'iter';
            options.TolFun = 1e-7;
            options.TolX = 1e-7;
            
            c_optimal = fminsearch(powerLawFit,c_initial_guesses,options);
            obj.power_law_fit_parameters = c_optimal;
            
            optimal_fit = getFullRadialPowerLawPred(c_optimal,size(data));
            obj.optimal_power_law_fit = optimal_fit;
            figure; subplot(1,3,1);
            nan_masked_data = data;
            nan_masked_data(data == -1) = nan;
            surf(nan_masked_data); shading flat;
            subplot(1,3,2);
            surf(log10(optimal_fit)); shading flat;
            subplot(1,3,3);
            residuals = data - optimal_fit;
            residuals(data == -1) = nan;
            surf(residuals); shading flat;
            colormap(flipud(bone));
        end
        
        
        
        
        function plotPowerLawMask(obj)
            % Plot the mask as a whole
            [data] = obj.singleLoad(2,2);
            exponent = 0.3;
            factor = 10;
            filtered_image = double(data).^exponent/factor;
            filtered_image(obj.power_law_fit_mask) = filtered_image(obj.power_law_fit_mask)*factor;
            % f2 = figure;
            % imagesc(filtered_image);
            f2 = plotDP(filtered_image,exponent);
            axis square;
            colormap jet
            colorbar
        end
        
        
        function fitDiskPositions(obj,use_BN_kernel)
            if ~use_BN_kernel
                probe = obj.loadProbe();
            else
                probe = obj.average_BN_probe;
            end
            stack = obj.partialLoad(1);
            
            stacksize = size(stack);
            probesize = size(probe);
            if stacksize(1) ~= probesize(1) || stacksize(2) ~= probesize(2)
                % assumes the scaling issue will be the same in each
                % direction.
                if use_BN_kernel  % pad with zeros
                    b1 = mean(probe(:,1));
                    b2 = mean(probe(:,end));
                    b3 = mean(probe(1,:));
                    b4 = mean(probe(end,:));
                    difference1 = stacksize(1) - probesize(1);
                    pad1num1 = ceil(difference1/2);
                    pad1num2 = floor(difference1/2);
                    padding1 = b1*ones(pad1num1,size(probe,2));
                    padding2 = b2*ones(pad1num2,size(probe,2));
                    probe = [padding1;probe;padding2];
                    difference2 = stacksize(2) - probesize(2);
                    pad2num1 = ceil(difference2/2);
                    pad2num2 = floor(difference2/2);
                    padding1 = b3*ones(size(probe,1),pad2num1);
                    padding2 = b4*ones(size(probe,1),pad2num2);
                    probe = [padding1,probe,padding2];
                else  % rescale a vacuum probe on bin 1
                    factor = stacksize(1)./probesize(1);
                    probe = imresize(probe,factor);
                end
                if size(probe,2) ~= stacksize(2)
                    probe(:,(stacksize(2)+1):end) = [];
                elseif size(probe,1) ~= stacksize(1)
                    probe(:,(stacksize(1)+1):end) = [];
                end
            end
            [sFit] = Gstrain01(stack,probe);
            sFit.posRefine = [];
%             [sFit] = Gstrain02(sFit);
            % Role of the function Gstrain02 has been replaced by 
            % setPositionFitLattices
            sFit.lat1 = [obj.beam_center_coords'; obj.lattice_vectors_G1];
            sFit.lat2 = [obj.beam_center_coords'; obj.lattice_vectors_G2];
            sFit.basis = [ones(size(obj.indices_G1,1),1) obj.indices_G1];
            sFit.xyInit1 = sFit.basis*sFit.lat1;
            sFit.xyInit2 = sFit.basis*sFit.lat2;
            sFit.xyFit1 = sFit.xyInit1;
            sFit.xyFit2 = sFit.xyInit2;
            
            for i = 1%:obj.num_load_chunks
                if i ~= 1
                    [stack,~,~] = obj.partialLoad(i);
                end
                sFit.stackSize = size(stack);
                [sFit] = Gstrain03(sFit);
                [sFit] = Gstrain04NPKmod(sFit,stack);
            end
            obj.sFit_struct = sFit;
            
% % %             sFit.stackSize(3) = obj.datacube_size(3);
            % Now these ones start doing plotting.
            [sFit] = Gstrain06(sFit);
            Gstrain07(sFit);
            Gstrain08(sFit);
            obj.sFit_struct = sFit;
            
        end
        
        
        
        % Originally copied and modified from setBlinkingDisks()
         function setPositionFitLattices(obj,explicit_chunks)
            if nargin < 2
                explicit_chunks = [];
            end
            
            if ~isempty(explicit_chunks)
                first_load = explicit_chunks(1);
                middle_load = explicit_chunks(2);
            else
                first_load = 1;
                middle_load = round(obj.num_load_chunks/2);
            end
            [data1,~,~] = obj.partialLoad(first_load);
            [data2,~,~] = obj.partialLoad(middle_load);
            diffraction_patterns = cat(3,data1,data2);
            
            [ disk_centers_G1, disk_centers_G2, indices_G1, indices_G2, lattice_vectors_G1, lattice_vectors_G2, beam_center_coords ] = ...
                makeDiskCentersForPositionFit( diffraction_patterns );
            obj.disk_centers_G1 = disk_centers_G1';
            obj.disk_centers_G2 = disk_centers_G2';
            obj.indices_G1 = indices_G1';
            obj.indices_G2 = indices_G2';
            obj.lattice_vectors_G1 = lattice_vectors_G1';
            obj.lattice_vectors_G2 = lattice_vectors_G2';
            obj.beam_center_coords = beam_center_coords';
            
         end
        
         
         
         % For when we don't have a vacuum probe to use.
         function [average_BN_probe,xcoords,ycoords] = setKernelFromHBN(obj,subtract_power_law)
             data = singleLoad(obj,2,2);
             plotDP(data);
             disp('Please zoom the graph until you are centered on hBN to define as kernel.');
             input('Press any key to continue.');
             disp('Please click twice on the graph to demarcate opposing corners of the square BN region.');
             [xcoords,ycoords] = ginput(2);
             xrange = ceil(min(xcoords):max(xcoords));
             yrange = ceil(min(ycoords):max(ycoords));
             % Let's assume that the BN doesn't shift and get the average
             % from that.
             aggregators = zeros(length(yrange),length(xrange),obj.num_load_chunks);
             for i = 1:obj.num_load_chunks
                 DPs = obj.partialLoad(i);
                 aggregators(:,:,i) = mean(mean(DPs(yrange,xrange,:,:),4),3);
             end
             obj.average_BN_probe = mean(mean(aggregators,4),3);
             average_BN_probe = obj.average_BN_probe;
         end
        
        
        
        
        %%%%%% Data manipulation methods %%%%%%
        
        % Chunk only in real space dimension 3.
        function [data,indices_start,indices_end] = partialLoad(obj,chunk_number)
            
            if chunk_number == obj.num_load_chunks  % Then need to be careful with the indices
                count = obj.datacube_size;
                m = mod(obj.datacube_size(3),obj.CHUNK_LOAD_WIDTH);
                if m == 0
                    count(3) = obj.CHUNK_LOAD_WIDTH;
                else
                    count(3) = m;
                end
            else
                count = obj.datacube_size;  % stays the same while start varies.
                count(3) = obj.CHUNK_LOAD_WIDTH;
            end
            %             fullsize = [1024,1024,40,40];
            fprintf('Loading chunk %d of %d.\n',chunk_number,obj.num_load_chunks);
            indices_start = [1,1,(chunk_number-1)*obj.CHUNK_LOAD_WIDTH+1,1];
            indices_end = indices_start + count - ones(1,4);
            data = h5read(obj.filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',indices_start,count);
            
        end
        
        
        
        function [data] = singleLoad(obj,i,j)
            start = [1,1,i,j];
            count = [obj.datacube_size(1),obj.datacube_size(2),1,1];
            data = h5read(obj.filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
        end
        
        
        
        function [probeimage,probeStruct] = loadProbe(obj)
            probeStruct = dm3Reader(obj.probe_filename);
            probeimage = probeStruct.image;
        end
        
        
        
        function plotProbe(obj)
            figure;
            probe = obj.loadProbe();
            pcolor(probe);
            colormap(gray);
            shading flat
        end
        
        
        
        function plotBNprobe(obj)
            figure;
            imagesc(obj.average_BN_probe);
            set(gca,'ydir','normal');
            colormap(gray);
            shading flat
            axis equal
        end
        
        
            
    end % End methods
    
    
    
    methods(Access = private)
        
        function darkFieldCallback(obj,position,plothandle)
            boolmat = obj.DSC_fit_storage(:,:,1) > position(1) & obj.DSC_fit_storage(:,:,1) < (position(1) + position(3)) ...
                & obj.DSC_fit_storage(:,:,2) > position(2) & obj.DSC_fit_storage(:,:,2) < (position(2) + position(4));
            axes(plothandle);
            %             boolmat = false(size(DSC1));
            h2 = imagesc(obj.xaxis,obj.yaxis,boolmat);
            colormap(gray);
            set(plothandle,'ydir','normal');
            axis equal
            xlim([obj.xaxis(1)-0.5*obj.scan_stepsize, obj.xaxis(end)+0.5*obj.scan_stepsize]);
            ylim([obj.yaxis(1)-0.5*obj.scan_stepsize, obj.yaxis(end)+0.5*obj.scan_stepsize]);
        end
        
    end
end

