classdef FourDSTEM_Analysis_Engine < handle
    % Class for performing the interferometry analysis of a 4DSTEM dataset
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
        disk_pixel_nums
        skipped_disks
        residuals_storage
        RMSR_storage
        DSC_fit_storage
        prefactor_fit_DSC_storage
        multistart_displacement_convergence  % for use when using
        classified_displacement_convergence
        scan_stepsize  % in units of nm, which is what is given by Maddie's filenames.
        exx  % strain map components
        eyy
        exy
        eyx
        gxy
        fixed_body_rotation
        FBR_type  % either 'radians' or 'degrees'
        theta_torsion
        vonMises
        NPK_strain_mag
        dfield_filtered_for_strain
        principal_strain_1
        principal_strain_2
        principal_angle
        principal_strain_max
        principal_strain_diff
        principal_strain_dilatation
        principal_strain_3Dcolor
        shear_axes
        xaxis
        yaxis
        xbase_for_strainmaps
        ybase_for_strainmaps
        rotated_basis
        initial_tensor_angle  % This is in degrees
        divide_by_two_for_intralayer
        sign_convention
        trim_value
        trim_tear_value
        elastic_fit_mask
        hBN_mask
        beamstop_mask
        graphene_mask
        detector_mask
        elastic_fit_parameters
        optimal_elastic_fit
        disk_centers_G1
        disk_centers_G2
        indices_G1
        indices_G2
        lattice_vectors_G1
        lattice_vectors_G2
        beam_center_coords  % This variable comes from position fit functions
        sFit_struct
        weight_vector
        average_BN_probe
        prefactor_storage
        trig_prefactors
        fitted_prefactors
        averaged_DP
        elastic_residuals_radial
        radial_correction_center
        elastic_residuals
        elastic_residuals_fit_parameters
        interpolated_residuals
        interpolated_polar_residuals
        strain_1D_roi_storage_cell
        strain_1D_percent_values
        soliton_width_1D_values
        AA_circlefit_values
        AA_circlefit_pixelerrors
        averageAA_DP
        soliton_walls_unique
        soliton_walls_merged
        SP_registration_unique_strainfiltermatched
        SP_registration_unique_dfield
        SP_registration_unique_filter
        untrimmed_AA_centers
        AB_area_storage  % sparse logicals, corresponding by number to the centroids below
        AB_centroids
        AB_segmentation_struct
        AB_label_matrix
        SP_transition_vectors
        emitter_pixel_pos
        emitter_displacements
        emitter_is_set
        annealed_dfield
        com_coords
        AA_unfiltered_weights
        AB_unfiltered_weights
        SP1_unfiltered_weights
        SP2_unfiltered_weights
        SP3_unfiltered_weights
        summed_weights
        ABangles
        ABsidelengths
        ABareas_nm
        ABcalc_boundary_edge_margin
        AAgaussianfit_FWHM_ellipses
        AAgaussianfit_raw_parameters
        AA_Gaussian_Circle_Fit_Params
        AAtriangulation
        triangle_sidelengths
        triangle_moire_angle
        moire_angle_estimates
        moire_angle_estimates_raw
        filter_immunity_mask
        uniaxial_moire_angles
        uniaxial_strain_angles
        uniaxial_strain_percents
        uniaxial_strain_residuals
        hotspot_suppression_mask
        tear_mask
        StrainStruct
        diffraction_theta_deg  % instance variable for rotational calibration (calculated)
        scan_direction_theta_deg  % instance variable for rotational calibration (external input)
    end
    
    methods
        % in the constructor, resultsfolderpath gives the full path to the
        % folder where images will be saved, including the name of the
        % folder itself.
        function obj = FourDSTEM_Analysis_Engine(filename,disk_filename,resultsfolderpath,scan_stepsize,probe_filename)
            % Assume that the filepath has been set in the driver
            % disk_filename is optional
            
            if nargin < 1
                return
            end
            
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
        
        
        
        
        
        % For use when the scan stepsize was not set correctly upon object
        % construction
        function reformAxes(obj)
            % Set x and y axes
            obj.xaxis = obj.scan_stepsize*(0:(obj.datacube_size(4)-1));
            obj.yaxis = obj.scan_stepsize*(0:(obj.datacube_size(3)-1));            
        end
        
        
        
        
        
        % Define the integration regions for the blinking calculations.
        function [] = setBlinkingDisks(obj,explicit_chunks,exponent)
            disp('Beginning the setBlinkingDisks() method...');
            if nargin < 2
                explicit_chunks = [];
            end
            if nargin < 3
                exponent = 0.2;
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
                [data1,~,~] = obj.partialLoad(first_load);
                [data2,~,~] = obj.partialLoad(middle_load);
                diffraction_patterns = cat(3,data1,data2);
                meanDP = mean(diffraction_patterns,3);
            elseif isempty(obj.averaged_DP)
                first_load = 1;
                middle_load = round(obj.num_load_chunks/2);
                [data1,~,~] = obj.partialLoad(first_load);
                [data2,~,~] = obj.partialLoad(middle_load);
                diffraction_patterns = cat(3,data1,data2);
                meanDP = mean(diffraction_patterns,3);
            else
                meanDP = obj.averaged_DP;
            end
            
            %             [ obj.disk_centers, obj.disk_radius ] = makeDiskCenters( diffraction_patterns );
            
            tf = input('Do you wish to mask off the hBN while selecting disks? 1/0');
            if tf
                this_hBN_mask = obj.hBN_mask;
            else
                this_hBN_mask = [];
            end
            %             if tf
            %                 disp('Prepare to mask off the hBN lattice and the beamstop.');
            %                 [ ~, this_hBN_mask, ~, beamstop_mask, ~ ] = makeHexagonalLatticeMask( meanDP );
            %
            %             else
            %                 this_hBN_mask = [];
            %             end
            
            while true
                [ disk_centers_temp, radius, skipped_disks ] = makeDiskCenters( meanDP, exponent, this_hBN_mask );
                % Removed by NPK on 03/15/2020 because we never use this.
                %                 randDPs_dim3 = randi(obj.CHUNK_LOAD_WIDTH,3,1);
                %                 randDPs_dim4 = randi(obj.datacube_size(4),3,1);
                %                 plotDP(diffraction_patterns(:,:,randDPs_dim3(1),randDPs_dim4(1)));
                %                 viscircles(disk_centers_temp,radius*ones(size(disk_centers_temp,1),1));
                %                 plotDP(diffraction_patterns(:,:,randDPs_dim3(2),randDPs_dim4(2)));
                %                 viscircles(disk_centers_temp,radius*ones(size(disk_centers_temp,1),1));
                %                 plotDP(diffraction_patterns(:,:,randDPs_dim3(3),randDPs_dim4(3)));
                %                 viscircles(disk_centers_temp,radius*ones(size(disk_centers_temp,1),1));
                
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
            obj.getNumBlinkingDiskPixels();
        end
        
        
        
        
        function getNumBlinkingDiskPixels(obj)
            xbase = 1:obj.datacube_size(1);
            ybase = 1:obj.datacube_size(2);
            [xspace,yspace] = meshgrid(xbase,ybase);
            obj.disk_pixel_nums = zeros(1,12);
            for q = 1:size(obj.disk_centers,1)
                x0 = obj.disk_centers(q,1);
                y0 = obj.disk_centers(q,2);
                tf = isInCircle(xspace,yspace,x0,y0,obj.disk_radius);
                obj.disk_pixel_nums(q) = nnz(tf);
            end
            obj.disk_pixel_nums(isnan(obj.disk_pixel_nums)) = 0;
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
        
        
        
        
        % NPK changed the way baselining was handled on 03/13/2020.
        % Subtracting off the radial residuals correction will now be
        % delegated to a helper function on a dataset-by-dataset manner,
        % and the baselined dataset will be integrated according to the
        % disk centers.
        function disk_average_storage = integrateDisks(obj,use_average_flag,subtractScattering,...
                subtractRadial,subtractCartesianInterpolation,...
                subtractPolarInterpolation,loadnum_crop_range)
            disp('Beginning the integrateDisks() method (sequential)...');
            if nargin < 7
                loadrange = 1:obj.num_load_chunks;
            else
                loadrange = loadnum_crop_range;
            end
            NUMDISKS = 12;
            disk_average_storage = zeros(obj.datacube_size(3),obj.datacube_size(4),NUMDISKS);
            center = obj.elastic_fit_parameters(1:2);
            
            for i = loadrange
                [thisdata,indices_start,indices_end] = obj.partialLoad(i);
                % Actually pull down the disk averages data that we wanted.
                [~,~,r,c] = size(thisdata);
                xbase = 1:obj.datacube_size(1);
                ybase = 1:obj.datacube_size(2);
                [xspace,yspace] = meshgrid(xbase,ybase);
                for rr = 1:r
                    for cc = 1:c
                        thisDP = thisdata(:,:,rr,cc);
                        if subtractScattering || subtractRadial || subtractCartesianInterpolation || subtractPolarInterpolation
                            % Delegate this slightly more involved calculation to the helper function.
                            %                             subtractScattering,subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation
                            %                             correctDP(obj,DP,subtract_elastic_flag,...
                            %                 radial_elastic_correction_flag,cartesian_interpolation_flag,polar_interpolation_flag)
                            thisDP = obj.correctDP(thisDP,subtractScattering,subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation);
                        end
                        %                             thisDP = double(thisDP) - obj.optimal_elastic_fit;
                        %                         elseif subtract_elastic_flag && radial_elastic_correction_flag
                        %
                        %
                        %                         end
                        global_row_index = indices_start(3) + rr - 1;
                        for q = 1:size(obj.disk_centers,1)
                            x0 = obj.disk_centers(q,1);
                            y0 = obj.disk_centers(q,2);
                            tf = isInCircle(xspace,yspace,x0,y0,obj.disk_radius);
                            %                             if radial_elastic_correction_flag
                            %                                 % compute radius based on indices
                            %                                 radius = sqrt((x0 - center(1)).^2 + (y0 - center(2)).^2);
                            %                                 [~,idx] = min(abs(obj.elastic_residuals_radial(:,1)-radius));
                            %                                 correction_val = obj.elastic_residuals_radial(idx,2);
                            %                             end
                            if use_average_flag
                                disk_average_storage(global_row_index,cc,q) = mean(thisDP(tf));
                            else
                                disk_average_storage(global_row_index,cc,q) = sum(thisDP(tf));
                            end
                        end
                    end
                end
                obj.disk_averages = disk_average_storage;
                clear thisdata
            end
        end
        
        
        function disk_average_storage = integrateDisksParallel(obj,subtract_power_law_flag)
            disp('Beginning the integrateDisks() method (parallel)...');
            NUMDISKS = 12;
            disk_average_storage = zeros(obj.datacube_size(3),obj.datacube_size(4),NUMDISKS);
            loopub = size(obj.disk_centers,1);
            plf = obj.optimal_elastic_fit;
            dcs = obj.disk_centers;
            dr = obj.disk_radius;
            for i = 1:obj.num_load_chunks
                [thisdata,indices_start,indices_end] = obj.partialLoad(i);
                % Actually pull down the disk averages data that we wanted.
                [~,~,r,c] = size(thisdata);
                xbase = 1:obj.datacube_size(1);
                ybase = 1:obj.datacube_size(2);
                [xspace,yspace] = meshgrid(xbase,ybase);
                local_disk_average_storage = zeros(r,c,size(obj.disk_centers,1));
                parfor rr = 1:r
                    for cc = 1:c
                        thisDP = thisdata(:,:,rr,cc);
                        if subtract_power_law_flag
                            thisDP = double(thisDP) - plf;
                        end
                        %                         global_row_index = indices_start(3) + rr - 1;
                        
                        for q = 1:loopub
                            tf = isInCircle(xspace,yspace,dcs(q,1),dcs(q,2),dr);
                            local_disk_average_storage(rr,cc,q) = mean(thisDP(tf));
                        end
                    end
                    
                    %                     if i == 1
                    %                         disk_average_storage = local_disk_average_storage;
                    %                     else
                    %                         disk_average_storage = cat(1,disk_average_storage,local_disk_average_storage);
                    %                     end
                end
                
                obj.disk_averages = disk_average_storage;
                clear thisdata
            end
        end
        
        
        % It looks like the right way to do this is to vectorize the
        % integration code. NPK
        function disk_average_storage = integrateDisksVectorized(obj,subtract_power_law_flag)
            disp('Beginning the integrateDisks() method (vectorized)...');
            NUMDISKS = 12;
            disk_average_storage = zeros(obj.datacube_size(3),obj.datacube_size(4),NUMDISKS);
            dr = obj.disk_radius;
            dc = obj.disk_centers;
            ds = obj.datacube_size;
            for i = 1:obj.num_load_chunks
                [thisdata,indices_start,indices_end] = obj.partialLoad(i);
                cs = indices_end(3)-indices_start(3)+1;
                
                % Actually pull down the disk averages data that we wanted.
                [~,~,r,c] = size(thisdata);
                xbase = 1:obj.datacube_size(1);
                ybase = 1:obj.datacube_size(2);
                [xspace,yspace] = meshgrid(xbase,ybase);
                %                 for rr = 1:r
                %                     for cc = 1:c
                %                         thisDP = thisdata(:,:,rr,cc);
                %                 if subtract_power_law_flag
                %                     thisDP = double(thisDP) - repmat(obj.optimal_elastic_fit,m4.datacube_size(1),m4.datacube_size(2),1,1);
                %                 end
                %                         global_row_index = indices_start(3) + rr - 1;
                for q = 1:size(obj.disk_centers,1)
                    x0 = dc(q,1);
                    y0 = dc(q,2);
                    tf = isInCircle(xspace,yspace,x0,y0,dr);
                    %                     disk_average_storage(global_row_index,cc,q) = mean(thisDP(tf));
                    toavpoints = thisdata(repmat(tf,1,1,cs,ds(4)));
                    yespts = nnz(tf);
                    r1 = reshape(toavpoints,[yespts,cs,ds(4)]);
                    r2 = mean(r1,1);
                    r3 = squeeze(r2);
                    disk_average_storage(indices_start(3):indices_end(3),indices_start(4):indices_end(4),q) = r3;
                end
            end
            %                 end
            obj.disk_averages = disk_average_storage;
            clear thisdata
            %             end
        end
        
        
        
        
        
        function figh = makeBlinkingPlots(obj,disk_ids,saveplotflag,shadingtype,normalize)
            if nargin < 2 || isempty(disk_ids)
                figh = makeMultiDataBlinkingPlot( obj.disk_averages(:,:,[1:3,7:9]),shadingtype,normalize,obj.xaxis,obj.yaxis);
            else
                figh = makeMultiDataBlinkingPlot( obj.disk_averages(:,:,disk_ids),shadingtype,normalize,obj.xaxis,obj.yaxis);
            end
            set(gcf, 'Position',  [0, 0, 1400, 800])
            if saveplotflag
                currentd = pwd;
                cd(obj.saved_plots_folderpath);
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
        function fitBlinking(obj,trig_prefactor,weight_vector_override,useParallel)
            
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
                weight_vector_local = ones(1,12);
            else
                weight_vector_local = weight_vector_override;
            end
            % Not removing the skipped disks will crash the program.
            weight_vector_local(obj.skipped_disks) = 0;
            weight_vector_local
            obj.weight_vector = weight_vector_local;
            
            %     weight_vector = [1,1,1,1,1,1,...
            %         0,0,0,0,0,0];
            %     weight_vector = [0,0,0,0,0,0,...
            %         1,1,1,1,1,1];
            
            %             makeMultiDataBlinkingPlot( disk_average_storage(:,:,[1:3,7:9]));
            nan_handle_flag = 0;
            for i = 1:obj.datacube_size(3)
                fprintf('Optimizing data row %d out of %d...\n',i,obj.datacube_size(3));
                if useParallel
                    parfor j = 1:obj.datacube_size(4)
                        this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
                        [DSCout,residOut,RMSRout] = obj.fitBlinkingInnards(this_disk_averages,trig_prefactor,nan_handle_flag,...
                            weight_vector_local,options);
                        DSC_fit_storage(i,j,:) = DSCout;
                        residuals_storage(i,j,:) = residOut;
                        RMSR_storage(i,j) = RMSRout;
                    end
                else
                    for j = 1:obj.datacube_size(4)
                        this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
                        [DSCout,residOut,RMSRout] = obj.fitBlinkingInnards(this_disk_averages,trig_prefactor,nan_handle_flag,...
                            weight_vector_local,options);
                        DSC_fit_storage(i,j,:) = DSCout;
                        residuals_storage(i,j,:) = residOut;
                        RMSR_storage(i,j) = RMSRout;
                    end
                end
            end
            
            obj.residuals_storage = residuals_storage;
            obj.RMSR_storage = RMSR_storage;
            obj.DSC_fit_storage = DSC_fit_storage;
        end
        
        
        
        
        % NPK 03/02/2020, fitting the trig prefactors globally across the
        % dataset.
        function fitBlinking2(obj)
            options = optimoptions('lsqnonlin');
            options.Display = 'iter';
            options.MaxFunctionEvaluations = 1000000;
            disk_average_storage = permute(obj.disk_averages,[1,3,2]);
            new_disk_average_storage = disk_average_storage(:,:,1);
            for i = 2:size(disk_average_storage,3)
                new_disk_average_storage = cat(1,new_disk_average_storage,disk_average_storage(:,:,i));
            end
            %             [ residuals ] = trigFittingFunction2outer( params, disk_average_storage )
            fitfun = @(params) trigFittingFunction2outer( params, new_disk_average_storage );
            diskmean = max(new_disk_average_storage);
            prefactor_relbounds = [0.5,1.5];
            
            lb = horzcat(diskmean*prefactor_relbounds(1),repmat([-1.24,0],1,size(new_disk_average_storage,1)));
            ub = horzcat(diskmean*prefactor_relbounds(2),repmat([1.24,1.43],1,size(new_disk_average_storage,1)));
            swarmsize = 500;
            for i = 1:swarmsize
                InitialPopulation(i,1).x = unifrnd(lb',ub');
            end
            Problem.ObjFunction = @(params) rms(fitfun(params'));
            Problem.Variables = numel(lb);
            Problem.LB = lb';
            Problem.UB = ub';
            opt = PSwarm('defaults');
            opt.MaxIter = 5E4;
            opt.MaxObj = 5E4;
            [BestParticle, BestParticleObj, RunData] = PSwarm(Problem,...
                InitialPopulation,opt);
            params_init = horzcat(mean(new_disk_average_storage),repmat([0,0.5],1,size(new_disk_average_storage,1)));
            %             results = lsqnonlin(fitfun,params_init,lb,ub,options);
            
            prefactor_results = results(1:12);
            DSC_x_results = results(13:2:(end-1));
            DSC_y_results = results(14:2:end);
            DSC_x_results = reshape(DSC_x_results,[obj.datacube_size(3),obj.datacube_size(4)]);
            DSC_y_results = reshape(DSC_y_results,[obj.datacube_size(3),obj.datacube_size(4)]);
            obj.DSC_fit_storage = cat(3,DSC_x_results,DSC_y_results);
            obj.prefactor_storage = prefactor_results;
        end
        
        
        
        
        
        % Realspace crop range is a cell array of two two-element vectors.
        function fitBlinking3(obj,weight_vector_override,useParallel,realspace_crop_range,useAnalyticJacobian,useFittedPrefactors)
            disp('Beginning the fitBlinking3() method...');
            if nargin < 2
                weight_vector_override = [];
                useParallel = 0;
                realspace_crop_range = [];
            end
            if (nargin < 5), useAnalyticJacobian = true; end
            if (nargin < 6 && ~isempty(obj.fitted_prefactors))
                useFittedPrefactors = true;
            elseif (nargin < 6 && ~isempty(obj.fitted_prefactors))
                useFittedPrefactors = false;
            end
            
            if isempty(obj.disk_averages)
                error('Please define the disks using obj.setBlinkingDisks()');
            end
            
            if ~isempty(realspace_crop_range)
                v1 = realspace_crop_range{1};
                v2 = realspace_crop_range{2};
                rng1 = v1(1):v1(2);
                rng2 = v2(1):v2(2);
                disk_average_storage = obj.disk_averages(rng1,rng2,:);
            else
                disk_average_storage = obj.disk_averages;
                rng1 = 1:obj.datacube_size(3);
                rng2 = 1:obj.datacube_size(4);
            end
            size1 = rng1(end)-rng1(1)+1;
            size2 = rng2(end)-rng2(1)+1;
            
            options = optimoptions('lsqnonlin');
            options.Display = 'off';
            if useAnalyticJacobian
                options.SpecifyObjectiveGradient = true;
                %                 options.CheckGradients = true;
            end
            
            DSC_fit_storage_local = zeros(size1,size2,2);
            residuals_storage_local = zeros(size1,size2,12);
            RMSR_storage_local = zeros(size1,size2);
            convergence_storage_local = cell(size1,size2);
            
            if isempty(weight_vector_override)
                weight_vector_local = ones(1,12);
            else
                weight_vector_local = weight_vector_override;
            end
            weight_vector_local(obj.skipped_disks) = 0;
            weight_vector_local
            obj.weight_vector = weight_vector_local;
            nan_handle_flag = 0;
            
            if ~useFittedPrefactors
                safe_trig_prefactors = obj.trig_prefactors;
            else
                safe_trig_prefactors = obj.fitted_prefactors;
            end
            safe_trig_prefactors(isnan(safe_trig_prefactors)) = 0;
            
            for i = 1:numel(rng1)
                fprintf('Optimizing data row %d out of %d...\n',i,size1);
                if useParallel
                    parfor j = 1:numel(rng2)
                        this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
                        [DSCout,residOut,RMSRout,convergence_storageOut] = obj.fitBlinkingInnards(this_disk_averages,safe_trig_prefactors,nan_handle_flag,...
                            weight_vector_local,options);
                        DSC_fit_storage_local(i,j,:) = DSCout;
                        residuals_storage_local(i,j,:) = residOut;
                        RMSR_storage_local(i,j) = RMSRout;
                        convergence_storage_local{i,j} = convergence_storageOut;
                    end
                else
                    % Changed 03/22 which should hopefully allow better
                    % processing of cropped data.
                    for j = 1:numel(rng2)
                        this_disk_averages = permute(disk_average_storage(i,j,:),[1,3,2]);
                        [DSCout,residOut,RMSRout,convergence_storageOut] = obj.fitBlinkingInnards(this_disk_averages,safe_trig_prefactors,nan_handle_flag,...
                            weight_vector_local,options);
                        DSC_fit_storage_local(i,j,:) = DSCout;
                        residuals_storage_local(i,j,:) = residOut;
                        RMSR_storage_local(i,j) = RMSRout;
                        convergence_storage_local{i,j} = convergence_storageOut;
                    end
                end
            end
            
            obj.residuals_storage = residuals_storage_local;
            obj.RMSR_storage = RMSR_storage_local;
            obj.DSC_fit_storage = DSC_fit_storage_local;
            obj.multistart_displacement_convergence = convergence_storage_local;
        end
        
        
        
        
        
        % Hamish-style fitting, specifically designed where there must be
        % an analytic Jacobian supplied to the fitting routine.
        %
        % Note that this fitting routine can only be run after a standard
        % fitting routine, because we need a good initial guess to try to
        % avoid a local minimum. However, the AA/AB discontinuous
        % interconversion may actually cause quite a bit of a problem here,
        % so keep an eye on it. Could even have separate fits for inner and
        % outer ring prefactor amplitudes. This could work quite well given
        % the contours of the inner residual surface -- no interconversion
        % barrier.
        function fitBlinking4(obj,weight_vector_override,useParallel,realspace_crop_range)
            % Assumption here of 12 disks being fit.
            if nargin < 4
                realspace_crop_range = [];
            end
            if ~isempty(realspace_crop_range)
                sz1 = realspace_crop_range(1);
                sz2 = realspace_crop_range(2);  % both of these will crop from 1 by default.
                DSC = obj.DSC_fit_storage(1:sz1,1:sz2,:);
                disk_avs = obj.disk_averages(1:sz1,1:sz2,:);
            else
                sz1 = obj.datacube_size(3);
                sz2 = obj.datacube_size(4);
                DSC = obj.DSC_fit_storage;
                disk_avs = obj.disk_averages;
                
            end
            if isempty(weight_vector_override)
                weight_vector_used = obj.weight_vector;
            else
                weight_vector_used = weight_vector_override;
            end
            assert(nnz(~(weight_vector_used == 1 | weight_vector_used == 0)) == 0,...
                'Full/sparse multiplication issues means only binary weight vectors can be used with fitBlinking4 at this time.');
            response_values_prep = reshape(disk_avs,[sz1*sz2,12])';
            response_values = response_values_prep(:);
            response_values(isnan(response_values)) = 0;
            objfcn = @(variables) HamishTrigFittingFunctionResidualsWrapper( variables, response_values, weight_vector_used );
            options = optimoptions('lsqnonlin');
            options.SpecifyObjectiveGradient = true;
            options.MaxIterations = 1000; % Testing run
            %             options.CheckGradients = true; % Testing run
            options.UseParallel = logical(useParallel);
            options.Display = 'iter';
            options.FunctionTolerance = 1e-5;
            options.StepTolerance = 1e-5;
            
            DSC_prep = reshape(DSC,[sz1*sz2,2])';
            DSC_prep = (DSC_prep(:))';
            init_guess = horzcat(obj.trig_prefactors,DSC_prep);
            init_RMSR = rms(objfcn(init_guess));
            
            results = lsqnonlin(objfcn,init_guess,[],[],options);
            
            final_RMSR = rms(objfcn(results));
            obj.fitted_prefactors = results(1:12);
            dscoutprep = reshape(results(13:end),[2,sz1*sz2]);
            dscoutx = dscoutprep(1,:);
            dscouty = dscoutprep(2,:);
            dscoutx = reshape(dscoutx,[sz1,sz2]);
            dscouty = reshape(dscouty,[sz1,sz2]);
            obj.prefactor_fit_DSC_storage = cat(3,dscoutx,dscouty);
            %
            %             figure
            %             getCustomDisplacementColor( displacement_to_use, [], [], 1 );
        end
        
        
        
        
        
        % Linear regression to obtain fractional character of AA, AB, SP1,
        % SP2, SP3
        function fitBlinking5(obj,use_fitted_prefactors,use_fnnls)
            %             [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
            if use_fitted_prefactors
                prefactors = obj.fitted_prefactors;
            else
                prefactors = obj.trig_prefactors;
            end
            AAbasis = ones(12,1);
            ABbasis = vertcat(0.25*ones(6,1),ones(6,1));
            SP1basis = zeros(12,1);
            SP1basis([3,6,8,11]) = 1;
            SP2basis = zeros(12,1);
            SP2basis([1,4,9,12]) = 1;
            SP3basis = zeros(12,1);
            SP3basis([2,5,7,10]) = 1;
            basis_vectors = prefactors'.*horzcat(AAbasis,ABbasis,SP1basis,SP2basis,SP3basis);
            basis_vectors(obj.skipped_disks,:) = [];
            [r,c,h] = size(obj.disk_averages);
            interferometry_values = reshape(obj.disk_averages,[r*c,h])';
            interferometry_values(obj.skipped_disks,:) = [];
            weights = zeros(size(basis_vectors,2),size(interferometry_values,2));
            if use_fnnls
                for j = 1:r*c
                    if mod(j,1000) == 0
                        fprintf('Row %d out of %d.\n',j,r*c);
                    end
                    y = interferometry_values(:,j);
                    XtX = basis_vectors'*basis_vectors;
                    Xty = basis_vectors'*y;
                    weights(:,j) = fnnls(XtX,Xty);
                end
            else
                weights = basis_vectors\interferometry_values;
            end
            raw_residuals = interferometry_values - basis_vectors*weights;
            AAweights = reshape(weights(1,:)',[r,c]);
            ABweights = reshape(weights(2,:)',[r,c]);
            SP1weights = reshape(weights(3,:)',[r,c]);
            SP2weights = reshape(weights(4,:)',[r,c]);
            SP3weights = reshape(weights(5,:)',[r,c]);
            summed_weights = AAweights + ABweights + SP1weights + SP2weights + SP3weights;
            rms_data = reshape(rms(interferometry_values)',[r,c]);
            percent_rms_residuals = reshape(rms(raw_residuals)',[r,c])./rms_data * 100;
            for i = 1:7
                switch i
                    case 1
                        plotmat = AAweights;
                    case 2
                        plotmat = ABweights;
                    case 3
                        plotmat = SP1weights;
                    case 4
                        plotmat = SP2weights;
                    case 5
                        plotmat = SP3weights;
                    case 6
                        plotmat = summed_weights;
                    case 7
                        plotmat = percent_rms_residuals;
                end
                figure;
                %                 if any(i == 1:5)
                set(gcf,'Position',[0,400,1100,400]);
                if i ~= 7
                    subplot(1,2,1);
                end
                imagesc(plotmat);
                colormap('fire');
                colorbar;
                axis equal
                set(gca,'yDir','normal');
                switch i
                    case 1
                        title('AA weights unnormalized');
                    case 2
                        title('AB weights unnormalized');
                    case 3
                        title('SP1 weights unnormalized');
                    case 4
                        title('SP2 weights unnormalized');
                    case 5
                        title('SP3 weights unnormalized');
                    case 6
                        title('Summed weights');
                    case 7
                        title('RMS residuals as % of data');
                end
                if i ~= 7
                    normalized_plotmat = plotmat./summed_weights * 100;
                    subplot(1,2,2);
                    imagesc(normalized_plotmat);
                    colormap('fire');
                    colorbar;
                    axis equal
                    set(gca,'yDir','normal');
                    switch i
                        case 1
                            title('AA weights % normalized');
                        case 2
                            title('AB weights % normalized');
                        case 3
                            title('SP1 weights % normalized');
                        case 4
                            title('SP2 weights % normalized');
                        case 5
                            title('SP3 weights % normalized');
                        case 6
                            title('Summed weights % normalized');
                    end
                end
            end
            obj.AA_unfiltered_weights = AAweights;
            obj.AB_unfiltered_weights = ABweights;
            obj.SP1_unfiltered_weights = SP1weights;
            obj.SP2_unfiltered_weights = SP2weights;
            obj.SP3_unfiltered_weights = SP3weights;
            obj.summed_weights = summed_weights;
        end
        
        
        
        
        
        % Provides a way to try to get rid of the AB hotspots directly
        function filtered_weights = linearDisplacementFilter1(obj,region,threshold,closeandfill,figureflag)
            if nargin < 3
                threshold = 0.3;
            end
            pixel_threshold = 40;
            if strcmp(region,'AA')
                filt = obj.AA_unfiltered_weights;
            elseif strcmp(region,'AB')
                filt = obj.AB_unfiltered_weights;
            else
                filt = region;
            end
            if figureflag
                figure;
                imagesc(filt);
                colormap('fire');
                colorbar;
                axis equal
                set(gca,'yDir','normal');
                filtered_weights = filt;
                title('Before area open');
            end
            
            logicals = filt > threshold;
            keep = bwareaopen(logicals,pixel_threshold);
            filt(~keep) = 0;
            
            if closeandfill
                
                filt = ~bwmorph(~filt,'clean',1);
                filt = ~bwmorph(~filt,'spur',3);
                filt = ~bwmorph(~filt,'hbreak',3);
                filt = ~bwneighborerode( ~filt, 5, 1, 0 );
                filt = ~bwmorph(~filt,'clean',1);
                filt = ~bwmorph(~filt,'spur',3);
                filt = ~bwmorph(~filt,'hbreak',3);
                filt = ~bwneighborerode( ~filt, 5, 1, 0 );
                filt = ~bwmorph(~filt,'clean',1);
                filt = ~bwmorph(~filt,'spur',3);
                filt = ~bwmorph(~filt,'hbreak',3);
                filt = ~bwneighborerode( ~filt, 5, 1, 0 );
                filt = ~bwmorph(~filt,'clean',1);
                filt = ~bwmorph(~filt,'spur',3);
                filt = ~bwmorph(~filt,'hbreak',3);
                filt = ~bwneighborerode( ~filt, 5, 1, 0 );
                %                 filt = bwmorph(filt,'thicken',5);
                %                 filt = ~bwneighborerode( ~filt, 5, 5, 0 );
            end
            
            filtered_weights = filt;
            
            if figureflag
                figure;
                imagesc(filt);
                colormap('fire');
                colorbar;
                axis equal
                set(gca,'yDir','normal');
                title('After area open');
            end
        end
        
        
        
        
        
        
        % Provides a way to identify AA regions, where hotspots should be
        % suppressed.
        function mask = linearDisplacementFilterMasks(obj,use_normalized,filtered_thresh,gausssigma,gaussthreshold,region,suppress_others,figureflag)
            if isa(region,'char')
                filtermat = obj.retrieveWeightArray(region,use_normalized);
            else
                filtermat = region;
            end
            
            if ~isempty(filtered_thresh)
                filtermat = obj.linearDisplacementFilter1(filtermat,filtered_thresh,0,figureflag);
            end
            
            % Cut others out if desired via recursion
            if suppress_others
                othermats = obj.retrieveWeightArray(sprintf('~%s',region),use_normalized);
                totalmask = false(size(filtermat));
                for i = 1:numel(othermats)
                    use_normalized = 1;
                    thismat = othermats{i};
                    thismask = obj.linearDisplacementFilterMasks(use_normalized,0.3,gausssigma,0.25,thismat,0,figureflag);
                    totalmask = totalmask | thismask;
                end
                AB_tailored = obj.linearDisplacementFilter1('AB',0.2,1,figureflag);
                totalmask = totalmask | AB_tailored;
                if figureflag
                    figure
                    imagesc(totalmask); axis equal; colormap(fire); colorbar; set(gca,'yDir','normal');
                end
                filtermat(totalmask) = 0;
            end
            
            filtered = imgaussfilt(filtermat,gausssigma);
            if figureflag
                figure;
                imagesc(filtered); axis equal; colormap(fire); colorbar; set(gca,'yDir','normal');
            end
            
            if suppress_others
                mask = false(size(filtered));
                fraction = 0.75;
                filtered(filtered < 0.01) = 0;
                domains = imregionalmax(filtered);
                domain_lininds = find(domains);
                MARGIN = 30;
                [xspace,yspace] = meshgrid(1:size(filtered,1),1:size(filtered,2));
                for i = 1:numel(domain_lininds)
                    maxheight = filtered(domain_lininds(i));
                    [domain_r,domain_c] = ind2sub(size(filtered),domain_lininds(i));
                    domainmask = true(size(filtered));
                    domainmask(xspace < domain_c-MARGIN | xspace > domain_c+MARGIN | yspace < domain_r-MARGIN | yspace > domain_r+MARGIN) = false;
                    mask = mask | (filtered > maxheight*fraction & domainmask);
                end
            else
                mask = filtered > gaussthreshold;
            end
            if figureflag
                figure; imagesc(mask); axis equal; set(gca,'yDir','normal');
            end
        end
        
        
        
        
        
        % These are good settings for dataset 24
        function mask = getHotspotSuppressionMask(obj)
            figureflag = false;
%             mask = obj.linearDisplacementFilterMasks(0,0.1,10,0.8,'AA',1,figureflag);
            mask = obj.linearDisplacementFilterMasks(1,0.1,10,0.8,'AA',1,figureflag);
        end
        
        
        
        
        
        function filtered_displacement = suppressHotspots(obj,displacement_to_use)
            METHOD = 1;
            x_disps = displacement_to_use(:,:,1);
            y_disps = displacement_to_use(:,:,2);
            mask = obj.getHotspotSuppressionMask();
            sub_lin_inds = find(mask);
            for q = sub_lin_inds'
                these_convergence_options = obj.classified_displacement_convergence{q};
                these_amplitude_options = sum(these_convergence_options{1}.^2,2).^0.5;
                [val,idx] = min(these_amplitude_options);
                if METHOD == 0
                    if val < 0.6
                        x_disps(q) = these_convergence_options{1}(idx,1);
                        x_disps(q) = these_convergence_options{1}(idx,2);
                    end
                elseif METHOD == 1
                    x_disps(q) = these_convergence_options{1}(idx,1);
                    x_disps(q) = these_convergence_options{1}(idx,2);
                end
                %                         y_disps(q) = newdisp(2);
                
                %                 numconvlocs = size(these_convergence_options{1},1);
                %                 for k = 2:numconvlocs  % clearly #1 didn't work out, which is why we find ourselves here
                %                     newdisp = these_convergence_options{1}(k,:);
                %                     new_disp_amp = sqrt(newdisp(1)^2 + newdisp(2)^2);
                %                     % Reevaluate the filter criterion
                %                     tf = abs(new_disp_amp - this_ampfilt) > threshold;
                %                     if ~tf  % then replace originals, we are under the threshold
                %                         x_disps(q) = newdisp(1);
                %                         y_disps(q) = newdisp(2);
                %                         break
                %                     end
                %                 end % If the criterion is never reached, we will never replace the displacement vector
            end
            filtered_displacement = cat(3,x_disps,y_disps);
        end
        
        
        
        
        
        
        % For very small twist-angle datasets such as 15 and 18, it's
        % important to not lose the AA region to the median filters while
        % still being able to clean up the rest of the dataset.
        function drawFilterImmunityRegions(obj)
            [f,ax] = obj.makeCustomDisplacementColorPlot([],[],[],[],1,0,0);
            disp('Begin drawing circular ROI regions which will not be affected by the built-in displacement filter.');
            count = 1;
            mask = false(obj.datacube_size(3:4));
            xbase = 1:obj.datacube_size(4);
            ybase = 1:obj.datacube_size(3);
            [xspace,yspace] = meshgrid(xbase,ybase);
            while true
                figure(f);
                input('Zoom to the desired region and press any key when ready.');
                fprintf('Click once for the center and second for the radius of the #%dth roi.\n',count);
                vals = ginput(2);
                r = sqrt(sum((vals(1,:)-vals(2,:)).^2));
                center = vals(1,:);
                [tf] = isInCircle(xspace,yspace,center(1),center(2),r);
                mask = mask | tf;
                count = count + 1;
                yn = input('Do you wish to draw another circular roi? 1/0');
                if ~yn
                    break
                end
            end
            obj.filter_immunity_mask = mask;
        end
        
        
        
        
        
        
        function outmat = retrieveWeightArray(obj,id,use_normalized)
            repo = {'AA','AB','SP1','SP2','SP3'};
            if id(1) == '~'
                tocheck = id(2:end);
                remove = strcmp(tocheck,repo);
                repo(remove) = [];
                for i = 1:numel(repo)
                    outmat{i} = obj.retrieveWeightArray(repo{i},use_normalized);
                end
                return
            end
            switch id
                case 'AA'
                    outmat = obj.AA_unfiltered_weights;
                case 'AB'
                    outmat = obj.AB_unfiltered_weights;
                case 'SP1'
                    outmat = obj.SP1_unfiltered_weights;
                case 'SP2'
                    outmat = obj.SP2_unfiltered_weights;
                case 'SP3'
                    outmat = obj.SP3_unfiltered_weights;
                case 'sum'
                    outmat = obj.summed_weights;
                    %                 case '~AA'
                    %                     outmat{1} = obj.AB_unfiltered_weights;
                    %                     outmat{2} = obj.SP1_unfiltered_weights;
                    %                     outmat{3} = obj.SP2_unfiltered_weights;
                    %                     outmat{4} = obj.SP3_unfiltered_weights;
            end
            if use_normalized
                outmat = outmat./obj.summed_weights;
            end
        end
        
        
        
        
        
        
        function classifyDisplacementConvergence(obj)
            xTol = 0.01;
            yTol = 0.01;
            fitTol = 0.01;
            [r,c] = size(obj.multistart_displacement_convergence);
            for i = 1:r
                fprintf('Classifying convergence for row %d of %d.\n',i,r);
                for j = 1:c
                    these_convloc_RMSRs = obj.multistart_displacement_convergence{i,j}{2};
                    these_convloc_displacements = obj.multistart_displacement_convergence{i,j}{1};
                    classified_RMSRs = [];
                    classified_displacements = [];
                    classified_counts = [];
                    % Evaluate each convergence location and bin it
                    % accordingly.
                    for q = 1:size(these_convloc_RMSRs,1)
                        this_RMSR = these_convloc_RMSRs(q);
                        this_displacement = these_convloc_displacements(q,:);
                        if isempty(classified_RMSRs) || isempty(classified_displacements)
                            classified_RMSRs(1,1) = these_convloc_RMSRs(q);
                            classified_displacements(1,1:2) = these_convloc_displacements(q,:);
                            classified_counts(1,1) = 1;
                        else
                            % iterate through existing convergence
                            % locations to see if there is a match
                            match_flag = 0;
                            for k = 1:size(classified_RMSRs,1)
                                classified_RMSR = classified_RMSRs(k);
                                classified_displacement = classified_displacements(k,:);
                                evalFittol = abs(classified_RMSR - this_RMSR) < fitTol;
                                evalDisptol = all(abs(classified_displacement - this_displacement) < fitTol);
                                if evalFittol && evalDisptol
                                    % update values if RMSR is better
                                    if this_RMSR < classified_RMSR
                                        classified_RMSRs(k,1) = this_RMSR;
                                        classified_displacements(k,:) = this_displacement;
                                    end
                                    classified_counts(k,1) = classified_counts(k,1) + 1;
                                    match_flag = 1;
                                    break
                                end
                            end
                            % if no match, then add a storage entry
                            if ~match_flag
                                classified_RMSRs(end+1,1) = this_RMSR;
                                classified_displacements(end+1,:) = this_displacement;
                                classified_counts(end+1,1) = 1;
                            end
                        end
                    end
                    augmat = [classified_RMSRs,classified_counts,classified_displacements];
                    prepmat = sortrows(augmat);
                    obj.classified_displacement_convergence{i,j} = {prepmat(:,3:4),prepmat(:,1),prepmat(:,2)};
                end
            end
        end
        
        
        
        
        function plotConvergenceFrequency(obj)
            counts = zeros(size(obj.classified_displacement_convergence));
            for i = 1:size(obj.classified_displacement_convergence,1)
                for j = 1:size(obj.classified_displacement_convergence,2)
                    counts(i,j) = obj.classified_displacement_convergence{i,j}{3}(1);
                end
            end
            totalnum = size(obj.multistart_displacement_convergence{1,1}{1},1);
            figure
            imagesc(counts);
            xlabel('X real space');
            ylabel('Y real space');
            axis equal
            colormap jet
            colorbar
            set(gca,'YDir','normal');
            title(sprintf('Convergence frequency of optimizer, %d multistarts used',totalnum));
            
            figure
            edges = 0.5:1:12.5;
            histogram(counts(:),edges);
            xlabel('# of optimizers converged to global minimum out of 12');
            ylabel('# of real space pixels');
            set(gca,'FontSize',14);
        end
        
        
        
        
        % Currently this is a median filter on the amplitude, replacing the
        % displacement with the next available value from the multistart
        % convergence list until either the filter threshold is met or the
        % list runs out.
        %
        % Modified by NPK on 05/19/2020 to allow the distinction between a
        % hard filter (replace no matter what) and a soft filter (replace
        % only if within the filter threshold). Also added inputs for
        % sequential usage in the new filtering convention.
        function filtered_displacement = displacementContinuityFilter(obj,threshold,range,hard_flag,xdisp_input,ydisp_input,use_hotspot_suppression)
            if nargin < 2
                threshold = [];
            end
            if isempty(threshold)
                threshold = 0;
            end
            if nargin < 3
                range = 3;
            end
            if nargin < 4
                hard_flag = false;
            end
            if nargin < 5
                x_disps = obj.DSC_fit_storage(:,:,1);
                y_disps = obj.DSC_fit_storage(:,:,2);
            else
                if ~isempty(xdisp_input)
                    x_disps = xdisp_input;
                    y_disps = ydisp_input;
                else
                    x_disps = obj.DSC_fit_storage(:,:,1);
                    y_disps = obj.DSC_fit_storage(:,:,2);
                end
            end
            if nargin < 7
                use_hotspot_suppression = false;
            end
            org_xdisps = x_disps;
            org_ydisps = y_disps;  % Saved in case of applying a filter immunity mask
                
            if isempty(range)
                range = 3;
            end
            if isempty(obj.classified_displacement_convergence)
                obj.classifyDisplacementConvergence();
            end
            
            amplitude = (x_disps.^2 + y_disps.^2).^0.5;
            ampfilt = medfilt2(amplitude,[range range],'symmetric');
            sub_amp = abs(amplitude - ampfilt) > threshold;
            sub_lin_inds = find(sub_amp);
            counts = 0;
            for q = sub_lin_inds'
                this_ampfilt = ampfilt(q);
                these_convergence_options = obj.classified_displacement_convergence{q};
                numconvlocs = size(these_convergence_options{1},1);
                criterion_values = zeros(numconvlocs-1,1);  % Added to implement the hard filter.
                for k = 1:numconvlocs  
                    % clearly #1 didn't work out, which is why we find ourselves here
                    % NPK note 05/19/2020: With the introduction of the
                    % hard filter, we have to compute all of these
                    % including the first to ensure that the first isn't
                    % discarded, even if pretty good.
                    newdisp = these_convergence_options{1}(k,:);
                    new_disp_amp = sqrt(newdisp(1)^2 + newdisp(2)^2);
                    % Reevaluate the filter criterion
                    criterion = abs(new_disp_amp - this_ampfilt);
                    criterion_values(k) = criterion;
                    if ~hard_flag && k ~= 1  % this shouldn't matter at all, because we know the threshold isn't met.
                        tf = criterion > threshold;
                        if ~tf  % then replace originals, we are under the threshold
                            x_disps(q) = newdisp(1);
                            y_disps(q) = newdisp(2);
                            counts = counts + 1;
                            break
                        end
                        % If the criterion is never reached, we will never replace the displacement vector
                    end  % Using hard filter; we will take whichever convergence location is closest to the median.
                    % So wait until all have been computed
                end
                if hard_flag
                    if isempty(criterion_values)  % indicating that there was only one convergence locations
                        continue
                    end
%                     if q == 27930
%                         disp('investigate');
%                     end
                    [~,k_idx] = min(criterion_values);
%                     k_idx = idx+1; % the indices are offset by one because we are going 2:numconvloc.
                    newdisp = these_convergence_options{1}(k_idx,:);
                    x_disps(q) = newdisp(1);
                    y_disps(q) = newdisp(2);
                    if k_idx ~= 1 % i.e. we chose any other convergence location
                        counts = counts + 1;
                    end
                end
            end
            
            if ~isempty(obj.filter_immunity_mask)
                x_disps(obj.filter_immunity_mask) = org_xdisps(obj.filter_immunity_mask);
                y_disps(obj.filter_immunity_mask) = org_ydisps(obj.filter_immunity_mask);
                warning('Filter immunity mask has been used!');
            end
            
            % 05/28/2020: this is going to draw from the
            % hotspot_suppression_mask instance variable, the analog of the
            % filter immunity mask.
            % Within this mask, we will set to the lowest amplitude
            % convergence location possible.
            if use_hotspot_suppression
                mask = obj.hotspot_suppression_mask;
                lininds = find(mask(:));
                for q = lininds'
                    these_convlocs = obj.classified_displacement_convergence{q};
                    these_displacements = these_convlocs{1};
                    these_amplitudes = sqrt(these_displacements(:,1).^2 + these_displacements(:,2).^2);
                    [~,idx] = min(these_amplitudes);
                    new_displacements = these_displacements(idx,:);
                    x_disps(q) = new_displacements(1);
                    y_disps(q) = new_displacements(2);
                end
                warning('Hotspot suppression mask has been used!');
            end
            
            filtered_displacement = cat(3,x_disps,y_disps);
            fprintf('displacementContinuityFilter() detected %d amplitude outliers at a threshold of %.2f and range of %.2f, making %d vector replacements.\n',nnz(sub_amp),threshold,range,counts);
        end
        
        
        
        
        % 05/28/2020: Interactive function for building a hotspot
        % suppression mask. Original purpose was for examining dataset 22,
        % which is very zoomed in but has low S/N and hotspots around the
        % edges.
        function setManualHotspotSuppressionMask(obj)
            amp = sqrt(sum(obj.DSC_fit_storage.^2,3));
            
            f = figure;
            imagesc(amp); set(gca,'yDir','normal'); axis equal; colormap gray;
            disp('Begin drawing polygon ROI regions which will be suppressed to minimum convergence amplitude by the built-in displacement filter.');
            count = 1;
            mask = false(obj.datacube_size(3:4));
%             xbase = 1:obj.datacube_size(4);
%             ybase = 1:obj.datacube_size(3);
%             [xspace,yspace] = meshgrid(xbase,ybase);
            while true
                figure(f);
                input('Zoom to the desired region and press any key when ready.');
                fprintf('Draw the polygon roi of the #%dth roi.\n',count);
                tf = roipoly();
                mask = mask | tf;
                count = count + 1;
                yn = input('Do you wish to draw another polygon roi? 1/0');
                if ~yn
                    break
                end
            end
            obj.hotspot_suppression_mask = mask;
        end
        
        
        
        
        
        
        %
        %
        %         function [f,ax1] = makeCustomDisplacementColorPlot(obj,SP_hsv,AB_hsv,continuity_threshold,filter_range)
        %         end
        %
        %         function [f,ax1] = makeCustomDisplacementColorPlotUniqueSP(obj,SP_hues,SP_saturations,continuity_threshold,filter_range)
        %         end
        
        % Instructions:
        %
        % If you want to colorize the different saddle points differently,
        % then set SP_unique_flag = 1, and pass a 3x2 matrix for both
        % SP_hsv and AB_hsv. Each row is a successive saddle point or AB
        % region. The hues should be the same between the
        function [f,ax1] = makeCustomDisplacementColorPlot(obj,SP_hsv,AB_hsv,continuity_threshold,filter_range,SP_unique_flag,usePrefactorFit,useHotspotSuppression,savePlot,useAnnealedDfield)
            if nargin < 2
                SP_hsv = [0.8,0.8,1];
                AB_hsv = [0.5,0.3,1];
            end
            if nargin < 4
                continuity_threshold = [];
            end
            if nargin < 5
                filter_range = 3;
            end
            if nargin < 6
                SP_unique_flag = 0;
            end
            if nargin < 7
                usePrefactorFit = 0;
            end
            if nargin < 8
                useHotspotSuppression = false;
            end
            if nargin < 9
                savePlot = false;
            end
            if nargin < 10
                useAnnealedDfield = false;
            end
            if isempty(SP_hsv) || isempty(AB_hsv)
                if ~SP_unique_flag
                    SP_hsv = [0.8,0.8,1];
                    AB_hsv = [0.5,0.3,1];
                else  % no coordination with the color triangle needed; let the function figure out its own colors.
                    SP_hsv = [];
                    AB_hsv = [];
                end
            end
            grid_density = 0.01;
            % [ raster_points ] = getDSCHexagonRaster(grid_density,2.461);
            xbase = -1.5:grid_density:1.5;
            ybase = 0:grid_density:1.5;
            [xspace,yspace] = meshgrid(xbase,ybase);
            catmat = cat(3,xspace,yspace);
            
            [ RGB_color_stack ] = getCustomDisplacementColor( catmat, SP_hsv, AB_hsv, SP_unique_flag );
            figure;
            imagesc(xbase,ybase,RGB_color_stack);   
            set(gca,'ydir','normal');
            axis equal
            title('2D Colormap');
            xlabel('x displacement');
            ylabel('y displacement');
            
            if isa(continuity_threshold,'struct')  % the new filter protocol
                if useAnnealedDfield
                    xdisp = obj.annealed_dfield(:,:,1);
                    ydisp = obj.annealed_dfield(:,:,2);
                else
                    xdisp = obj.DSC_fit_storage(:,:,1);
                    ydisp = obj.DSC_fit_storage(:,:,2);
                end
                filterstruct = continuity_threshold;
                plotflag = 0;
                [ xdisp_filt,ydisp_filt ] = filterDisplacement( xdisp,ydisp,filterstruct,plotflag,obj );
                displacement_to_use = cat(3,xdisp_filt,ydisp_filt);
            else
                if usePrefactorFit
                    displacement_to_use = obj.prefactor_fit_DSC_storage;
                else
                    if ~isempty(continuity_threshold)
                        if isempty(obj.hotspot_suppression_mask) || ~useHotspotSuppression
                            displacement_to_use = obj.displacementContinuityFilter(continuity_threshold,filter_range);
                        else
                            hard_flag = 1;  % by default here.
                            displacement_to_use = obj.displacementContinuityFilter(continuity_threshold,filter_range,hard_flag,[],[],useHotspotSuppression);
                        end
                    else
                        displacement_to_use = obj.DSC_fit_storage;
                    end
                    if useHotspotSuppression && isempty(obj.hotspot_suppression_mask)
                        displacement_to_use = obj.suppressHotspots(displacement_to_use);
                    end
                end
            end
            
            if useAnnealedDfield  % Need to convert back to reduced zone
                xdisp = displacement_to_use(:,:,1);
                ydisp = displacement_to_use(:,:,2);
                [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( [xdisp(:),ydisp(:)] );
                xdisp_rz = reshape(reduced_zone_disps(:,1),size(xdisp));
                ydisp_rz = reshape(reduced_zone_disps(:,2),size(ydisp));
                displacement_to_use = cat(3,xdisp_rz,ydisp_rz);
            end
            
            [ RGB_color_stack ] = getCustomDisplacementColor( displacement_to_use, SP_hsv, AB_hsv, SP_unique_flag );
            f = figure;
            set(f,'Position',[200,200,800,700]);
            if ~SP_unique_flag
                ax1 = subplot(1,2,1);
                ax2 = subplot(1,2,2);
                axes(ax1);
            else
                ax1 = gca;
            end
            % NPK changed on 05/24/2020 from the following to make dataset 10 work:
%             if numel(obj.xaxis) == size(displacement_to_use,1) && numel(obj.yaxis) == size(displacement_to_use,2)
            if numel(obj.xaxis) == size(displacement_to_use,2) && numel(obj.yaxis) == size(displacement_to_use,1)
                imagesc(obj.xaxis,obj.yaxis,RGB_color_stack);
                xlabel('Real space x (nm)','FontSize',12);
                ylabel('Real space y (nm)','FontSize',12);
            else
                imagesc(RGB_color_stack);
                xlabel('Pixels x','FontSize',12);
                ylabel('Pixels y','FontSize',12);
                warning('Data dimensions do not seem to match axis dimensions. This may be due to fitting a cropped portion of the full datacube.');
            end
            set(ax1,'ydir','normal');
            axis square
            
            title('Colorized AB, SP, and AA regions','FontSize',14);
            hold on
            if ~SP_unique_flag
                set(ax1,'Position',[0.05 0.1 0.7 0.75]);
            end
            set(ax1,'FontSize',12);
            
            if ~SP_unique_flag
                [ RGB_color_stack_legend ] = getTriangleColorLegend(grid_density,AB_hsv,SP_hsv);
                legaxpos = [0.8,0.1,0.15,0.15];
                set(ax2,'Position',legaxpos);
                axes(ax2);
                imagesc(RGB_color_stack_legend);
                set(ax2,'ydir','normal');
                set(ax2,'Box','off');
                set(ax2,'XMinorTick','off');
                set(ax2,'YMinorTick','off');
                set(ax2,'XTick',[]);
                set(ax2,'YTick',[]);
                set(ax2,'visible','off');
                
                vertical_increment = 0.06;
                horizontal_increment = 0;% 0.0125;
                axsize = [0.05,0.05];
                th1 = axes('Position',[legaxpos(1)-horizontal_increment, legaxpos(2)-vertical_increment/1.4,axsize]);
                text(.025,.6,'AA','FontSize',12,'HorizontalAlignment','Center')
                set(th1,'visible','off');
                th1.XTick = [];
                th1.YTick = [];
                th2 = axes('Position',[legaxpos(1)+legaxpos(3)-horizontal_increment, legaxpos(2)-vertical_increment/1.4,axsize]);
                text(.025,.6,'SP','FontSize',12,'HorizontalAlignment','Center')
                set(th2,'visible','off');
                th2.XTick = [];
                th2.YTick = [];
                th3 = axes('Position',[legaxpos(1)+legaxpos(3)/2-horizontal_increment,legaxpos(2)+legaxpos(4)-vertical_increment/2,axsize]);
                text(.025,.6,'AB','FontSize',12,'HorizontalAlignment','Center')
                set(th3,'visible','off');
                th3.XTick = [];
                th3.YTick = [];
                th4 = axes('Position',[legaxpos(1)+legaxpos(3)/2-horizontal_increment,legaxpos(2)+legaxpos(4)+vertical_increment/4,axsize]);
                set(th4,'visible','off');
                th4.XTick = [];
                th4.YTick = [];
                text(.025,.6,'Stacking Order:','FontSize',14,'HorizontalAlignment','Center')
            end
            
            axes(ax1);
            f = f;
            
            if savePlot
                currentd = pwd;
                cd(obj.saved_plots_folderpath);
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                savefig(f,'2DColormapDisplacementField');
                saveas(f,'2DColormapDisplacementField.png');
                cd(currentd);
            end
        end
        
        
        
        
        
        % Function for using the custom 2D colormap from outside of the
        % class. Useful when using the filterDisplacement() /
        % filterStrain() revised function series.
        %
        % The assumption will be that the desired filtering operations have
        % already been completed.
        %
        % Nathanael Kazmierczak, 05/19/2020
        function [figh,axh] = makeOutsideCustomDisplacementColorPlot(obj,xdisp,ydisp)
            old_DSC = obj.DSC_fit_storage;
            obj.DSC_fit_storage = cat(3,xdisp,ydisp); % hacky but works
            [figh,axh] = obj.makeCustomDisplacementColorPlot([],[],[],[],1,0,0);
            obj.DSC_fit_storage = [];
            obj.DSC_fit_storage = old_DSC;
        end
        
        
        
        
        
        
        % These used to be in with fitBlinking() itself.
        %
        % Modified 07/18/2020: median_filter_threshold can now handle a
        % filterstruct passed in instead of a two-component vector, which
        % will get routed to filterDisplacement.m.
        function [fighs] = makeDisplacementMapsFromBlinking(obj,saveplotflag,median_filter_threshold,individual_residual_flag,shading_type)
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
                if isstruct(median_filter_threshold)
                    filterstruct = median_filter_threshold;
                    [ DSC1, DSC2 ] = filterDisplacement( DSC1,DSC2,filterstruct,false,obj );
                else
                    [ DSC1, DSC2 ] = npkMedFilt2Displacement( DSC1, DSC2, median_filter_threshold(1), median_filter_threshold(2) );
                end
            end
            scatter(DSC1(:),DSC2(:),3,'Filled');
            axh = gca;
            try
                plotFullDisplacementHexagons(axh);
            catch
                warning('Failed to plot full displacement hexagons in obj.makeDisplacementMapsFromBlinking()');
            end
            xlim([-1.5,1.5]);
            ylim([-0.1,1.5]);
            xlabel('Cartesian displacement component 1');
            ylabel('Cartesian displacement component 2');
            title('Displacement vector scatterplot');
            
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
            %             colormap(redblue)
            colormap(fire);
            axis equal
            shading interp
            colorbar
            title(sprintf(strcat(['Displacement field amplitude (angstroms)'])));
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
            title(sprintf(strcat(['Displacement field angle (radians)'])));
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
            
            figh_DSCresiduals_imagesc = figure;
            imagesc(RMSR_storage);
            shading 'flat'
            colormap(jet);
            colorbar
            xlabel('x');
            ylabel('y');
            title('Fit RMSR');
            
            if individual_residual_flag
                diskres_fighs = zeros(1,12);
                disk_strings = cell(1,12);
                disk_strings_imagesc = cell(1,12);
                disk_strings_figs = cell(1,12);
                disk_strings_figs_imagesc = cell(1,12);
                diskres_imagesc_fighs = zeros(1,12);
                for i = 1:12
                    if obj.weight_vector(i)
                        figh = figure;
                        this_disk_residual = obj.residuals_storage(:,:,i);
                        surfc(obj.xaxis,obj.yaxis,this_disk_residual);
                        title(sprintf('Interferometry Residuals for Disk %d',i));
                        xlabel('Realspace nm');
                        ylabel('Realspace nm');
                        shading flat
                        colormap(jet);
                        colorbar;
                        diskres_fighs(i) = figh;
                        disk_strings{i} = sprintf('Disk%dResidualsPlot.png',i);
                        disk_strings_figs{i} = sprintf('Disk%dResidualsPlot',i);
                        figh2 = figure;
                        imagesc(obj.xaxis,obj.yaxis,this_disk_residual);
                        title(sprintf('Interferometry Residuals for Disk %d',i));
                        xlabel('Realspace nm');
                        ylabel('Realspace nm');
                        colormap(jet);
                        colorbar
                        diskres_imagesc_fighs(i) = figh2;
                        disk_strings_imagesc{i} = sprintf('Disk%dResidualsPlotImageSC.png',i);
                        disk_strings_figs_imagesc{i} = sprintf('Disk%dResidualsPlotImageSC',i);
                    end
                end
            end
            
            fighs = [figh_DSCscatter,figh_DSCamp,figh_DSCangle,figh_DSCresiduals];
            
            if saveplotflag
                currentd = pwd;
                cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                saveas(figh_DSCscatter,'DSCscatter.png');
                saveas(figh_DSCamp,'DSCamplitude.png');
                saveas(figh_DSCangle,'DSCangle.png');
                saveas(figh_DSCresiduals,'DSCresiduals.png');
                saveas(figh_DSCresiduals_imagesc,'DSCresidualscontour.png');
                savefig(figh_DSCscatter,'DSCscatter');
                savefig(figh_DSCamp,'DSCamplitude');
                savefig(figh_DSCangle,'DSCangle');
                savefig(figh_DSCresiduals,'DSCresiduals');
                savefig(figh_DSCresiduals_imagesc,'DSCresidualscontour');
                if individual_residual_flag
                    for i = 1:12
                        if obj.weight_vector(i)
                            saveas(diskres_fighs(i),disk_strings{i});
                            savefig(diskres_fighs(i),disk_strings_figs{i});
                            saveas(diskres_imagesc_fighs(i),disk_strings_imagesc{i});
                            savefig(diskres_imagesc_fighs(i),disk_strings_figs_imagesc{i});
                        end
                    end
                end
                cd(currentd);
            end
            
        end
        
        
        
        
        
        % Prepare alternate coordinate systems for strain mapping
        % use the annealed displacement field by default because that's all
        % we care about from a strain mapping perspective.
        function [transformed_adisps,transformed_coords,xh1,xh2,yh1,yh2] = changeDisplacementBasis(obj)
            % (1) Prepare coordinate basis in terms of the Moiré
            xbase = (((1:obj.datacube_size(3))-1)/obj.scan_stepsize)+1;
            ybase = (((1:obj.datacube_size(4))-1)/obj.scan_stepsize)+1;
            %             ybase = 1:obj.datacube_size(4);
            [xspace,yspace] = meshgrid(xbase,ybase);
            for i = 1:2
                this_soliton = obj.soliton_walls_unique{i};
                soliton_seg = bwconncomp(full(this_soliton));
                num_walls = soliton_seg.NumObjects;
                xcoords_all = [];
                ycoords_all = [];
                for j = 1:num_walls
                    this_wall = soliton_seg.PixelIdxList{j};
                    xcoords_centered = xspace(this_wall) - mean(xspace(this_wall));
                    ycoords_centered = yspace(this_wall) - mean(yspace(this_wall));
                    % These will be linear so we can get the angle
                    % directly. Could also do a least squares fit but this
                    % might be fine.
                    xcoords_all = [xcoords_all;xcoords_centered];
                    ycoords_all = [ycoords_all;ycoords_centered];
                end
                m = xcoords_all\ycoords_all;
                % use the unit vector
                s(i,:) = [1/sqrt(1+m^2),m/sqrt(1+m^2)];
            end
            s1 = s(1,:);
            s2 = s(2,:);
            [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
            v1 = v1/norm(v1);
            v2 = v2/norm(v2);
            S = [s1',s2'];
            V = [v1',v2'];
            % for testing
            V = [1/sqrt(2), -1/sqrt(2);
                1/sqrt(2), 1/sqrt(2)]; % should be an orthonormal basis.
            
            adisps1 = obj.annealed_dfield(:,:,1);
            adisps1 = adisps1(:);
            adisps2 = obj.annealed_dfield(:,:,2);
            adisps2 = adisps2(:);
            adisps = vertcat(adisps1',adisps2');
            transformed_adisps_lin = (V*adisps)';
            transformed_adisps = cat(3,reshape(transformed_adisps_lin(:,1),size(xspace)),reshape(transformed_adisps_lin(:,2),size(xspace)));
            xcoords = xspace(:);
            ycoords = yspace(:);
            Cmat = vertcat(xcoords',ycoords');
            transformed_coords_lin = (S*Cmat)';
            transformed_coords = cat(3,reshape(transformed_coords_lin(:,1),size(xspace)),reshape(transformed_coords_lin(:,2),size(xspace)));
            
            xh1 = transformed_coords(2,1,1) - transformed_coords(1,1,1);
            xh2 = transformed_coords(1,2,1) - transformed_coords(1,1,1);
            yh1 = transformed_coords(2,1,2) - transformed_coords(1,1,2);
            yh2 = transformed_coords(1,2,2) - transformed_coords(1,1,2);
        end
        
        
        
        
        
         % Prepare alternate coordinate systems for strain mapping
        % use the annealed displacement field by default because that's all
        % we care about from a strain mapping perspective.
        %
        % Modifications made by NPK on 06/23/2020 to get a more reliable
        % angle quantification. Small pieces of soliton wall are taken out
        % and the edges of the image are cropped.
        
        function [transformed_coords,rotated_basis,xh1,xh2,yh1,yh2,angle,basisrotangle] = changeDisplacementBasis2(obj,SPnum_for_xaxis,trimvalue,minsolitonpixelsize)
            % (1) Prepare coordinate basis in terms of the Moiré
            obj.reformAxes();
            xbase = (((1:obj.datacube_size(4))-1)/obj.scan_stepsize)+1;  
            ybase = (((1:obj.datacube_size(3))-1)/obj.scan_stepsize)+1;
            % changed from the following on 05/24/2020 to make dataset 10
            % work.
%             xbase = (((1:obj.datacube_size(3))-1)/obj.scan_stepsize)+1;  
%             ybase = (((1:obj.datacube_size(4))-1)/obj.scan_stepsize)+1;
            if nargin < 3  % then take off 1/10th of the image from each side
                nside = numel(xbase);
                trimvalue = round(nside/10);  % This is in pixels, not nm.
            end
            if nargin < 4
                minsolitonpixelsize = 10;  
                % Heuristic value for how long a soliton needs to be before it provides a reliable contribution to the angle information 
            end
            
            %             ybase = 1:obj.datacube_size(4);
            [xspace,yspace] = meshgrid(xbase,ybase);
            if nargin < 2
                % Obtain the registration of the purple soliton (the
                % original method)
                this_soliton = obj.soliton_walls_unique{1};
            else
                % Added 06/03/2020, this new functionality allows the user
                % to choose which soliton wall direction will be used for
                % setting the x-axis.
                this_soliton = obj.soliton_walls_unique{SPnum_for_xaxis};
            end
            % Trim edges of soliton and filter small components if desired
            soliton_to_use = full(this_soliton);
            if ~isempty(trimvalue)
                soliton_to_use = trimArray(soliton_to_use,trimvalue,false(size(soliton_to_use)));
            end
            if ~isempty(minsolitonpixelsize)
                soliton_to_use = bwareaopen(soliton_to_use,minsolitonpixelsize);
            end
            soliton_seg = bwconncomp(soliton_to_use);
            num_walls = soliton_seg.NumObjects;
            xcoords_all = [];
            ycoords_all = [];
            for j = 1:num_walls
                this_wall = soliton_seg.PixelIdxList{j};
                xcoords_centered = xspace(this_wall) - mean(xspace(this_wall));
                ycoords_centered = yspace(this_wall) - mean(yspace(this_wall));
                % These will be linear so we can get the angle
                % directly. Could also do a least squares fit but this
                % might be fine.
                xcoords_all = [xcoords_all;xcoords_centered];
                ycoords_all = [ycoords_all;ycoords_centered];
            end
            m = xcoords_all\ycoords_all;
            % use the unit vector
            s = [1/sqrt(1+m^2),m/sqrt(1+m^2)];
            angle = atan(s(2)/s(1));  % The extra rotation needs to make this -pi/2
            newrotangle = -pi/2 - angle;
            newrotangle = -newrotangle;  % Because we are rotating the coordinate basis in the other direction, rather than rotating the saddle point
            basisrotangle = newrotangle;
            coords_rotmat = [cos(newrotangle), -sin(newrotangle); sin(newrotangle), cos(newrotangle)];
            coords_rotmat2 = [cos(-newrotangle), -sin(-newrotangle); sin(-newrotangle), cos(-newrotangle)];
            
            
%             s1 = s(1,:);
%             s2 = s(2,:);
%             [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
%             v1 = v1/norm(v1);
%             v2 = v2/norm(v2);
%             S = [s1',s2'];
% % %             V = [v1',v2'];
% % %             % for testing
% % %             V = [1/sqrt(2), -1/sqrt(2);
% % %                 1/sqrt(2), 1/sqrt(2)]; % should be an orthonormal basis.
% % %             
% % %             adisps1 = obj.annealed_dfield(:,:,1);
% % %             adisps1 = adisps1(:);
% % %             adisps2 = obj.annealed_dfield(:,:,2);
% % %             adisps2 = adisps2(:);
% % %             adisps = vertcat(adisps1',adisps2');
% % %             transformed_adisps_lin = (V*adisps)';
% % %             transformed_adisps = cat(3,reshape(transformed_adisps_lin(:,1),size(xspace)),reshape(transformed_adisps_lin(:,2),size(xspace)));
            xcoords = xspace(:);
            ycoords = yspace(:);
            Cmat = vertcat(xcoords',ycoords');
            transformed_coords_lin = (coords_rotmat2*Cmat)';
            transformed_coords = cat(3,reshape(transformed_coords_lin(:,1),size(xspace)),reshape(transformed_coords_lin(:,2),size(xspace)));
            
            basis = [1 0; 0 1];
            rotated_basis = coords_rotmat*basis;
            
            xh1 = transformed_coords(2,1,1) - transformed_coords(1,1,1);
            xh2 = transformed_coords(1,2,1) - transformed_coords(1,1,1);
            yh1 = transformed_coords(2,1,2) - transformed_coords(1,1,2);
            yh2 = transformed_coords(1,2,2) - transformed_coords(1,1,2);
            
            if false
            figure; imagesc(xspace); set(gca,'yDir','normal'); axis equal; colorbar; title('Original x coordinates');
            hold on
            quiver([100,100],[150,150],basis(1,:)*10,basis(2,:)*10);
            
            figure; imagesc(transformed_coords(:,:,1)); set(gca,'yDir','normal'); axis equal; colorbar; title('Transformed x coordinates');
            hold on
            quiver(100,100,s(1)*10,s(2)*10);
            quiver([100,100],[150,150],rotated_basis(1,:)*10,rotated_basis(2,:)*10);
            end
        end
        
        
        
        
        
        % Obtaining the numerical gradient.
        function makeStrainMapsFromBlinking(obj,saveplotflag,shadingtype,colormaptype,trim_edge_pixel_num,median_filter_threshold,TVpre_fraction,TVpost_fraction,filterstruct)
            if nargin < 3
                median_filter_threshold = 0;
            end
            if isempty(obj.DSC_fit_storage)
                error('Run fitBlinking() before makeStrainMaps().');
            end
            NUM_CONTOUR_LINES = 50;
            % DSC_fit_storage is a 3-way array in which the first component
            % of the 3rd dimension is ux(x,y) and the second is uy(x,y)
            %             uX = obj.DSC_fit_storage(:,:,1);
            %             uY = obj.DSC_fit_storage(:,:,2);
            % 04/03/2020 NPK changed this to the annealed vector field,
            % because otherwise strain maps are meaningless.
            
            uX = obj.annealed_dfield(:,:,1);
            uY = obj.annealed_dfield(:,:,2);
            
            
            [ uX,uY ] = filterDisplacement( uX,uY,filterstruct,2,obj );
            
            changeOfBasis = 0;
            METHOD = 'old';
            MEAN = true;
            FILTER = false;
%             MEDFILTTHRES = true;
            SMOOTH = false;
            scale_color_map_flag = true;
            vonMises = true;
            PAD = false;
            ZEROFILL = false;
            
            if changeOfBasis > 0
                [transformed_adisps,transformed_coords,xh1,xh2,yh1,yh2] = obj.changeDisplacementBasis();
                uX = transformed_adisps(:,:,1);
                uY = transformed_adisps(:,:,2);
            end
            
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            if trim_edge_pixel_num(1) > 0
                xbase((end-trim_edge_pixel_num(1)+1):end) = [];
                xbase(1:trim_edge_pixel_num(1)) = [];
                ybase((end-trim_edge_pixel_num(1)+1):end) = [];
                ybase(1:trim_edge_pixel_num(1)) = [];
                uX = trimArray(uX,trim_edge_pixel_num(1));
                uY = trimArray(uY,trim_edge_pixel_num(1));
            end
            
            if PAD
%                 uX = padarray(uX,[50,50],'replicate','both');
%                 uY = padarray(uY,[50,50],'replicate','both');
            end
            
            if FILTER
                orguX = uX;
                orguY = uY;
                uX = medfilt2(uX);
                uY = medfilt2(uY);
                % Sets the corners to zero, annoyingly, so we must correct
                % that
                uX = trimArray(uX,1,orguX);
                uY = trimArray(uY,1,orguY);
            end
%             if MEDFILTTHRES
%             end
            
            plotflag = 1;
            if TVpre_fraction(1) > 0 || TVpre_fraction(2) > 0
                uX = TVM_image_denoise( uX, TVpre_fraction(1), plotflag, 'X displacement', colormaptype );
                uY = TVM_image_denoise( uY, TVpre_fraction(2), plotflag, 'Y displacement', colormaptype );
            end
            
            if ~isempty(median_filter_threshold)
                [ uX, uY ] = npkMedFilt2Displacement( uX, uY, median_filter_threshold(1), median_filter_threshold(2) );
            end
            
            
            
            spacing = obj.scan_stepsize;  % This is in nm and we want to keep it that way.
            %             uX = 0.1*uX;  % converts the X displacement to nm
            %             uY = 0.1*uY;  % converts the X displacement to nm
            touse = cat(3,uX,uY);
            
            
            if strcmp(METHOD,'old')
                if changeOfBasis == 1
                    [obj.exx,obj.exy] = gradient(uX,spacing);
                    [obj.eyx,obj.eyy] = gradient(uY,spacing);
                elseif changeOfBasis == 2
                    [obj.exx,obj.exy] = gradient(uX,xh1,xh2);
                    [obj.eyx,obj.eyy] = gradient(uY,yh1,yh2);
                else
                    [obj.exx,obj.exy] = gradient(uX,spacing);
                    [obj.eyx,obj.eyy] = gradient(uY,spacing);
                end
            else
                %                 touse(:,:,1) = imgaussfilt(medfilt2(touse(:,:,1)),0.001);
                %                 touse(:,:,2) = imgaussfilt(medfilt2(touse(:,:,2)),0.001);
                permutation_number = 5;
                PAD = 2;
                xbounds = zeros(1,2);
                ybounds = zeros(1,2);
                xbounds(1) = max(obj.emitter_displacements(:,1)) + PAD;
                xbounds(2) = min(obj.emitter_displacements(:,1)) - PAD;
                ybounds(1) = max(obj.emitter_displacements(:,2)) + PAD;
                ybounds(2) = min(obj.emitter_displacements(:,2)) - PAD;
                [obj.exx,obj.exy] = computeGradient(touse,'x',permutation_number,xbounds,ybounds);
                [obj.eyx,obj.eyy] = computeGradient(touse,'y',permutation_number,xbounds,ybounds);
                %                 obj.exx = imgaussfilt(obj.exx,5);
            end
            
            % Convert into units of nm/nm *100 (strain %) for the displacement
            % rather than angstroms/nm (what it is calculated in, and the
            % only thing that computeGradient will take)
            obj.exx = 10*obj.exx;
            obj.exy = 10*obj.exy;
            obj.eyx = 10*obj.eyx;
            obj.eyy = 10*obj.eyy;
            
            
            if FILTER
%                 obj.exx = medfilt2(medfilt2(obj.exx));
%                 obj.exy = medfilt2(medfilt2(obj.exy));
%                 obj.eyx = medfilt2(medfilt2(obj.eyx));
%                 obj.eyy = medfilt2(medfilt2(obj.eyy));
                exxorg = obj.exx;
                exyorg = obj.exy;
                eyxorg = obj.eyx;
                eyyorg = obj.eyy;
                obj.exx = medfilt2(obj.exx);
                obj.exy = medfilt2(obj.exy);
                obj.eyx = medfilt2(obj.eyx);
                obj.eyy = medfilt2(obj.eyy);
                obj.exx  = trimArray(obj.exx,1,exxorg);
                obj.exy  = trimArray(obj.exy,1,exyorg);
                obj.eyx  = trimArray(obj.eyx,1,eyxorg);
                obj.eyy  = trimArray(obj.eyy,1,eyyorg);
            end
            
            
            if trim_edge_pixel_num(2) > 0
                xbase((end-trim_edge_pixel_num(2)+1):end) = [];
                xbase(1:trim_edge_pixel_num(2)) = [];
                ybase((end-trim_edge_pixel_num(2)+1):end) = [];
                ybase(1:trim_edge_pixel_num(2)) = [];
                obj.exx = trimArray(obj.exx,trim_edge_pixel_num(2));
                obj.exy = trimArray(obj.exy,trim_edge_pixel_num(2));
                obj.eyx = trimArray(obj.eyx,trim_edge_pixel_num(2));
                obj.eyy = trimArray(obj.eyy,trim_edge_pixel_num(2));
            end
            
            
            if MEAN
                obj.exx = obj.exx - mean(mean(obj.exx));
                obj.eyy = obj.eyy - mean(mean(obj.eyy));
                obj.exy = obj.exy - mean(mean(obj.exy));
                obj.eyx = obj.eyx - mean(mean(obj.eyx));
            end
            
            if ZEROFILL
                if PAD
                    obj.exx = padarray(obj.exx,[0,50],'symmetric','both');
                    obj.exy = padarray(obj.exy,[0,50],'symmetric','both');
                    obj.eyx = padarray(obj.eyx,[0,50],'symmetric','both');
                    obj.eyy = padarray(obj.eyy,[0,50],'symmetric','both');
                    obj.exx = padarray(obj.exx,[50,0],'symmetric','both');
                    obj.exy = padarray(obj.exy,[50,0],'symmetric','both');
                    obj.eyx = padarray(obj.eyx,[50,0],'symmetric','both');
                    obj.eyy = padarray(obj.eyy,[50,0],'symmetric','both');
%                     addones = ones(size(obj.exx));
%                     obj.exx = trimArray(obj.exx,50,mean(mean(obj.exx))*addones);
%                     obj.exy = trimArray(obj.exy,50,mean(mean(obj.exy))*addones);
%                     obj.eyx = trimArray(obj.eyx,50,mean(mean(obj.eyx))*addones);
%                     obj.eyy = trimArray(obj.eyy,50,mean(mean(obj.eyy))*addones);
                else
                    error('look at me');
                end
            end
            
            % TVM filter which actually works
            if TVpost_fraction > 0
                obj.exx = TVM_image_denoise( obj.exx, TVpost_fraction, plotflag, 'exx strain', colormaptype );
                obj.exy = TVM_image_denoise( obj.exy, TVpost_fraction, plotflag, 'exy strain', colormaptype );
                obj.eyx = TVM_image_denoise( obj.eyx, TVpost_fraction, plotflag, 'eyx strain', colormaptype );
                obj.eyy = TVM_image_denoise( obj.eyy, TVpost_fraction, plotflag, 'eyy strain', colormaptype );
            end
            
            
            
            
            if SMOOTH
                sval = 0.5;
                obj.exx = imgaussfilt(obj.exx,sval);
                obj.exy = imgaussfilt(obj.exy,sval);
                obj.eyx = imgaussfilt(obj.eyx,sval);
                obj.eyy = imgaussfilt(obj.eyy,sval);
            end
            % NPK 04/27/2020 block used to be here before zero filling
            % experiments
%             if MEAN
%                 obj.exx = obj.exx - mean(mean(obj.exx));
%                 obj.eyy = obj.eyy - mean(mean(obj.eyy));
%                 obj.exy = obj.exy - mean(mean(obj.exy));
%                 obj.eyx = obj.eyx - mean(mean(obj.eyx));
%             end
            obj.gxy = obj.exy + obj.eyx;
            obj.theta_torsion = 0.5*(obj.eyx - obj.exy);
            
            if strcmp(METHOD,'new')
                xbase = 1:80;
                ybase = 1:80;
            end
            
            if trim_edge_pixel_num(3) > 0
                xbase((end-trim_edge_pixel_num(3)+1):end) = [];
                xbase(1:trim_edge_pixel_num(3)) = [];
                ybase((end-trim_edge_pixel_num(3)+1):end) = [];
                ybase(1:trim_edge_pixel_num(3)) = [];
                obj.exx = trimArray(obj.exx,trim_edge_pixel_num(3));
                obj.exy = trimArray(obj.exy,trim_edge_pixel_num(3));
                obj.eyx = trimArray(obj.eyx,trim_edge_pixel_num(3));
                obj.eyy = trimArray(obj.eyy,trim_edge_pixel_num(3));
                obj.gxy = trimArray(obj.gxy,trim_edge_pixel_num(3));
                obj.theta_torsion = trimArray(obj.theta_torsion,trim_edge_pixel_num(3));
            end
            
            
            % Using the following function to trim edges (which may have
            % artifacts)
            %             [ A_trimmed ] = trimArray( A,trim_width )
            figh_exx = figure;
            plot_exx = obj.exx; % trimArray(obj.exx,trim_edge_pixel_num);
            
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_exx);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_exx,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('exx strain % (nm x displacment / nm x real space * 100)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            figh_eyy = figure;
            plot_eyy = obj.eyy;%trimArray(obj.eyy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_eyy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_eyy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('eyy strain % (nm y displacment / nm y real space * 100)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            figh_gxy = figure;
            plot_gxy = obj.gxy;%trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_gxy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_gxy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('gxy total shear strain % (nm x(y) displacment / nm y(x) real space * 100)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            figh_theta = figure;
            plot_theta = obj.theta_torsion;%trimArray(obj.theta_torsion,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_theta);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_theta,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            colorbar
            title('Torsional strain %: theta (nm x(y) displacment / nm y(x) real space * 100)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            
            if vonMises
                vM = obj.exy.*obj.eyx - obj.exx.*obj.eyy;
                figure;
                imagesc(xbase,ybase,vM);
                %                 scaleColorMap( colormaptype, 0);
                colormap(colormaptype);
                axis equal
                xlim([xbase(1) xbase(end)]);
                ylim([ybase(1) ybase(end)]);
                colorbar
                title('vonMises Strain Invariant');
                xlabel('x (nm)');
                ylabel('y (nm)');
                set(gca,'ydir','normal');
            end
            
            
            
            
            if saveplotflag
                currentd = pwd;
                cd(obj.saved_plots_folderpath)
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
        
        
        
        
        
        
        
        
         % Obtaining the numerical gradient.
         % 05/20/2020: by default, this will rotate the basis so that the
         % purple saddle point is facing straight up and down.
         %
        % sign_convention argument is for subtracting the estimated moire
        % angle from the shear strain maps. It appears that different
        % annealing methods may make these come out differently (e.g. -/-
        % means on DS 26 and 4, but +/+ means on DS 7), so we will do it
        % manually for each dataset. It should be pretty obvious what the
        % correct convention is in each case -- could even be automated.
        function [eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                makeStrainMapsFromBlinking2(obj,saveplotflag,colormaptype,trim_edge_pixel_num,...
            filterstruct,sign_convention,trim_tear_pixel_num,overlay_registration,divide_by_two_for_intralayer,SP_num_for_xaxis)
            
            if (nargin < 7), trim_tear_pixel_num = []; end;
            if (nargin < 8), overlay_registration = false; end
            if (nargin < 9), divide_by_two_for_intralayer = true; end;
            if (nargin < 10), SP_num_for_xaxis = 1; end;  % corresponds to the purple saddle points we have been consistently examining.
            
            NUM_CONTOUR_LINES = 50;
            ALPHA_CONSTANT = 0.25;
            shadingtype = 'flat';
            scale_color_map_flag = true;
%             vonMises = true;
            
            % 04/03/2020 NPK changed this to the annealed vector field,
            % because otherwise strain maps are meaningless.
            uX = obj.annealed_dfield(:,:,1);
            uY = obj.annealed_dfield(:,:,2);
            
            [ uX,uY ] = filterDisplacement( uX,uY,filterstruct,0,obj );
            [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( [uX(:),uY(:)] );
            uXr = trimArray(reshape(reduced_zone_disps(:,1),size(uX)),trim_edge_pixel_num(1));
            uYr = trimArray(reshape(reduced_zone_disps(:,2),size(uX)),trim_edge_pixel_num(1));
            
            if ~isempty(obj.tear_mask) && ~isempty(trim_tear_pixel_num)
                if ~isempty(obj.tear_mask)
                    tear_boundary_mask = boundarymask(obj.tear_mask);
                    tear_boundary_mask_thick = tear_boundary_mask;
                    for i = 1:trim_tear_pixel_num
                        tear_boundary_mask_thick = boundarymask(tear_boundary_mask_thick) | tear_boundary_mask_thick;
                        tear_boundary_mask_thick(obj.tear_mask) = false;
                    end
                end
                uX(tear_boundary_mask_thick) = nan;
                uY(tear_boundary_mask_thick) = nan;
            end
            
            [figh_filtadisp,~] = makeOutsideCustomDisplacementColorPlot(obj,uXr,uYr);
            title('Reduced zone representation of filtered data for strain mapping');
            % NPK 06/03/2020: This has been changed so that there is an
            % option for using an SP value other than one for setting the
            % x-axis.
            if SP_num_for_xaxis ~= 0
                [transformed_coords,rotated_basis,xh1,xh2,yh1,yh2] = changeDisplacementBasis2(obj,SP_num_for_xaxis);
            else
                % Principally useful for testing, to make sure unrotated
                % strain maps match those from Colin. Directional
                % derivative work will come later.
                rotated_basis = eye(2);
            end
            
%             obj.makeOutsideCustomDisplacementColorPlot(uX,uY);
%             % Add the basis plot
            
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            if trim_edge_pixel_num(1) > 0
                xbase((end-trim_edge_pixel_num(1)+1):end) = [];
                xbase(1:trim_edge_pixel_num(1)) = [];
                ybase((end-trim_edge_pixel_num(1)+1):end) = [];
                ybase(1:trim_edge_pixel_num(1)) = [];
                uX = trimArray(uX,trim_edge_pixel_num(1));
                uY = trimArray(uY,trim_edge_pixel_num(1));
            end
            
            
            % Rework the units here, so it is nm on the top and the bottom
            % from the beginning. Otherwise, application of the derived
            % intralayer strain formulas can't be done.
            spacing = obj.scan_stepsize;  % This is in nm and we want to keep it that way.
            uX = 0.1*uX;  % converts the X displacement to nm
            uY = 0.1*uY;  % converts the y displacement to nm
            
            % Compute the intralayer strain from the interlayer
            % displacement field. This section uses results from the page
            % "Final Strain Formulas" in Nathanael's binder, 05/22/2020
            changeOfBasis = 3;
            if changeOfBasis == 3
                % First, get the gradient in the original x,y system of the
                % pixels.
                [dUinterX_dorgx,dUinterX_dorgy] = gradient(uX,spacing);
                [dUinterY_dorgx,dUinterY_dorgy] = gradient(uY,spacing);
                % Next, compute the directional derivatives of the x and y
                % components of the interlayer strain. Do this by dot
                % products with a unit normal vector, using the directional
                % derivative formula.
                xdir = rotated_basis(:,1)';
                ydir = rotated_basis(:,2)';
                dUinterX = vertcat(dUinterX_dorgx(:)',dUinterX_dorgy(:)');
                dUinterX_dx = reshape((xdir*dUinterX)',size(uX));
                dUinterX_dy = reshape((ydir*dUinterX)',size(uX));
                dUinterY = vertcat(dUinterY_dorgx(:)',dUinterY_dorgy(:)');
                dUinterY_dx = reshape((xdir*dUinterY)',size(uX));
                dUinterY_dy = reshape((ydir*dUinterY)',size(uX));
                % Now, these derivatives represent the interlayer strain in
                % the new coordinate system. But we want to convert them to
                % intralayer strain using the moire twist angle, estimated
                % independently.
                if ~isempty(obj.moire_angle_estimates)
                    mangles = trimArray(obj.moire_angle_estimates,trim_edge_pixel_num(1));  % Will need to crop this to match.
                else
                    warning('There is no Moire angle set. This could lead to inaccurate shear strain maps.');
                    tf = input('Do you wish to proceed? 1/0:');
                    if ~tf
                        error('Aborting program');
                    end
                    mangles = zeros(size(dUinterX_dx));
                end
                
                % Old, prior to 05/31/2020 modification.
%                 obj.exx = dUinterX_dx;
%                 obj.exy = dUinterX_dy + sign_convention(1)*deg2rad(mangles);
%                 obj.eyx = dUinterY_dx + sign_convention(2)*deg2rad(mangles);
%                 obj.eyy = dUinterY_dy;
                
                % 05/31/2020: realized that the division by two was
                % missing, so these were net strain maps, if you will. This
                % affects all previous, and makes sense (i.e., we were
                % getting numbers near 0.8% for shear soliton strain just
                % as Muller calculated, but each layer should only bear
                % 0.4% of that strain).
                %
                % NPK 06032020 DOUBLE CHECK THIS. I thought Muller's
                % derivation already accounted for the fact that there are
                % two graphene layers, and divided by two.
                
                
                if divide_by_two_for_intralayer
                    obj.exx = 0.5*(dUinterX_dx);
                    obj.exy = 0.5*(dUinterX_dy + sign_convention(1)*deg2rad(mangles));
                    obj.eyx = 0.5*(dUinterY_dx + sign_convention(2)*deg2rad(mangles));
                    obj.eyy = 0.5*(dUinterY_dy);
                else
                    obj.exx = (dUinterX_dx);
                    obj.exy = (dUinterX_dy + sign_convention(1)*deg2rad(mangles));
                    obj.eyx = (dUinterY_dx + sign_convention(2)*deg2rad(mangles));
                    obj.eyy = (dUinterY_dy);
                end
                
                % 05/30/2020: reworked handling of the cross derivative
                % terms after figuring out the issue with the vonMises
                % strain invariant.
                %
                % Note that these values are just reparameterizations of
                % the simple shear strain. They are not to be multiplied by
                % 1/2 again; that is simply part of the definition for the
                % strain tensor.
                obj.gxy = 0.5*(obj.exy + obj.eyx);
                obj.fixed_body_rotation = 0.5*(obj.eyx - obj.exy);
            end
            
            
            
            % Convert into units of nm/nm *100 (strain %) for the displacement
            % rather than nm/nm (what it is computed in after the revision)
            % Note 05/30/2020. Be VERY careful with using the cross derivatives directly.
            % Fixed body rotation is convoluted into them. It may prove
            % that we need to just use pure shear strain when talking about
            % soliton walls? But then again, following Muller's definition,
            % it may not.
            obj.exx = 100*obj.exx;
            obj.exy = 100*obj.exy;
            obj.eyx = 100*obj.eyx;
            obj.eyy = 100*obj.eyy;
            obj.gxy = 100*obj.gxy;
            % 05/30/2020: express fixed-body rotation in terms of degrees.
            obj.fixed_body_rotation = rad2deg(obj.fixed_body_rotation);
           
            
            % Install more sophisticated filters here
            
            
            if trim_edge_pixel_num(2) > 0
                xbase((end-trim_edge_pixel_num(2)+1):end) = [];
                xbase(1:trim_edge_pixel_num(2)) = [];
                ybase((end-trim_edge_pixel_num(2)+1):end) = [];
                ybase(1:trim_edge_pixel_num(2)) = [];
                obj.exx = trimArray(obj.exx,trim_edge_pixel_num(2));
                obj.exy = trimArray(obj.exy,trim_edge_pixel_num(2));
                obj.eyx = trimArray(obj.eyx,trim_edge_pixel_num(2));
                obj.eyy = trimArray(obj.eyy,trim_edge_pixel_num(2));
            end
            
% OLD, removed 05/30/2020            
%             obj.gxy = obj.exy + obj.eyx;
%             obj.theta_torsion = 0.5*(obj.eyx - obj.exy);
            % The new total strain metric: sqrt of sum of squares. Still
            % has units of nm/nm %
            % Modified 05/30/2020 to remove fixed-body rotation.
%             obj.NPK_strain_mag = sqrt(obj.exx.^2 + obj.exy.^2 + obj.eyy.^2 + obj.eyx.^2); 
            obj.NPK_strain_mag = sqrt(obj.exx.^2 + obj.gxy.^2 + obj.eyy.^2); 
                        
            if trim_edge_pixel_num(3) > 0
                xbase((end-trim_edge_pixel_num(3)+1):end) = [];
                xbase(1:trim_edge_pixel_num(3)) = [];
                ybase((end-trim_edge_pixel_num(3)+1):end) = [];
                ybase(1:trim_edge_pixel_num(3)) = [];
                obj.exx = trimArray(obj.exx,trim_edge_pixel_num(3));
                obj.exy = trimArray(obj.exy,trim_edge_pixel_num(3));
                obj.eyx = trimArray(obj.eyx,trim_edge_pixel_num(3));
                obj.eyy = trimArray(obj.eyy,trim_edge_pixel_num(3));
                obj.gxy = trimArray(obj.gxy,trim_edge_pixel_num(3));
                obj.theta_torsion = trimArray(obj.theta_torsion,trim_edge_pixel_num(3));
            end
            
            
            figh_DSCscatter = figure;
            DSC1 = uX;
            DSC2 = uY;
            scatter(DSC1(:)*10,DSC2(:)*10,3,'Filled');
            axh = gca;
            plotFullDisplacementHexagons(axh);
%             xlim([-1.5,1.5]);
%             ylim([-0.1,1.5]);
            xlabel('Cartesian displacement component 1');
            ylabel('Cartesian displacement component 2');
            title('Displacement vector scatterplot after pre-strain filter');
            
            
            xcolor = 'g';
            ycolor = 'm';
            
            
            %%%%% exx
            figh_exx = figure;
            plot_exx = obj.exx; % trimArray(obj.exx,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_exx);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_exx,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('exx intralayer normal strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('exx interlayer normal strain');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            
            %%%%% eyy
            figh_eyy = figure;
            plot_eyy = obj.eyy;%trimArray(obj.eyy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_eyy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_eyy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('eyy intralayer normal strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('eyy interlayer normal strain');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            %%%%% exy figure
            figh_exy = figure;
            plot_exy = obj.exy; %trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_exy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_exy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('exy intralayer simple shear strain (cross derivative)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('exy interlayer simple shear strain (cross derivative)');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            %%%%% eyx figure
            figh_eyx = figure;
            plot_eyx = obj.eyx;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_eyx);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_eyx,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('eyx intralayer simple shear strain (cross derivative)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('eyx interlayer simple shear strain (cross derivative)');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            
            %%%%% vonMises strain
            % 05/30/2020: the revised definition of vonMises, having
            % removed the fixed-body rotation.
%             vM = obj.exy.*obj.eyx - obj.exx.*obj.eyy;
            vM = obj.gxy.^2 - obj.exx.*obj.eyy;
            figh_vonMises = figure;
            imagesc(xbase,ybase,vM);
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('Intralayer vonMises Strain Invariant (fixed body rotation removed)');
                ylabel(cbh,'(Intralayer strain %)^2');
            else
                title('Interlayer vonMises Strain Invariant (fixed body rotation removed)');
                ylabel(cbh,'(Interlayer strain %)^2');
            end
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            
            
            
            %%%% gxy pure shear strain plot
            figh_gxy = figure;
            plot_gxy = obj.gxy;%trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_gxy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_gxy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('gyx intralayer pure shear strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('gyx interlayer pure shear strain');
                ylabel(cbh,'Interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            
            %%%% Fixed body rotation figure
            figh_theta = figure;
            plot_theta = obj.fixed_body_rotation;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_theta);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_theta,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if sign_convention(1) == 0 && sign_convention(2) == 0
                if divide_by_two_for_intralayer
                    title('Intralayer fixed body rotation, Moire rotation included');
                    ylabel(cbh,'Intralayer total rotation (degrees)');
                else
                    title('Interlayer fixed body rotation, Moire rotation included');
                    ylabel(cbh,'Interlayer total rotation (degrees)');
                end
            else
                if divide_by_two_for_intralayer
                    title('Intralayer fixed body rotation, Moire rotation removed');
                    ylabel(cbh,'Intralayer reconstruction rotation (degrees)');
                else
                    title('Interlayer fixed body rotation, Moire rotation removed');
                    ylabel(cbh,'Interlayer reconstruction rotation (degrees)');
                end
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            %%%% NPK total strain figure.  % Note: should this include gxy
            % twice?
            figh_gradcomponentmag = figure;
            plot_mag = obj.NPK_strain_mag;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_mag);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_mag,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            colormap(fire);  % a diverging colormap makes less sense here since it will be all positive
            cbh = colorbar;
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            if divide_by_two_for_intralayer
                title('Root-mean-square intralayer strain magnitude of exx, eyy, and gxy');
                ylabel(cbh,'Strain %');
            else
                title('Root-mean-square interlayer strain magnitude of exx, eyy, and gxy');
                ylabel(cbh,'Strain %');
            end
%             title('Total strain metric (3-vector magnitude of all strain components but rotation) (nm/nm %)');
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
            
            %%%% Make statistics
            eyx_percentiles = prctile(obj.eyx(:),0:100);
            gradientsum_percentiles = prctile(obj.NPK_strain_mag(:),0:100);
            eyx_values = obj.eyx;
            gradientsum_values = obj.NPK_strain_mag;
            
            if saveplotflag
                currentd = pwd;
                cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                
                if overlay_registration
                    saveas(figh_filtadisp,'FilteredDispForGradient_overlay.png');
                    saveas(figh_exx,'exx_strain_overlay.png');
                    saveas(figh_eyy,'eyy_strain_overlay.png');
                    saveas(figh_exy,'exy_strain_overlay.png');
                    saveas(figh_eyx,'eyx_strain_overlay.png');
                    saveas(figh_gxy,'gxy_pure_shear_strain_overlay.png');
                    saveas(figh_theta,'fixed_body_rotation_overlay.png');
                    saveas(figh_vonMises,'vonMises_strain_overlay.png');
                    saveas(figh_gradcomponentmag,'GradientComponentsMagnitude_overlay.png');
                    saveas(figh_DSCscatter,'StrainFilterDisplacementScatterplot.png');
                    
                    savefig(figh_filtadisp,'FilteredDispForGradient_overlay');
                    savefig(figh_exx,'exx_strain_overlay');
                    savefig(figh_eyy,'eyy_strain_overlay');
                    savefig(figh_exy,'exy_strain_overlay');
                    savefig(figh_eyx,'eyx_strain_overlay');
                    savefig(figh_gxy,'gxy_pure_shear_strain_overlay');
                    savefig(figh_theta,'fixed_body_rotation_overlay');
                    savefig(figh_vonMises,'vonMises_strain_overlay');
                    savefig(figh_gradcomponentmag,'GradientComponentsMagnitude_overlay');
                    savefig(figh_DSCscatter,'StrainFilterDisplacementScatterplot');
                else
                    saveas(figh_filtadisp,'FilteredDispForGradient.png');
                    saveas(figh_exx,'exx_strain.png');
                    saveas(figh_eyy,'eyy_strain.png');
                    saveas(figh_exy,'exy_strain.png');
                    saveas(figh_eyx,'eyx_strain.png');
                    saveas(figh_gxy,'gxy_pure_shear_strain.png');
                    saveas(figh_theta,'fixed_body_rotation.png');
                    saveas(figh_vonMises,'vonMises_strain.png');
                    saveas(figh_gradcomponentmag,'GradientComponentsMagnitude.png');
                    saveas(figh_DSCscatter,'StrainFilterDisplacementScatterplot.png');
                    
                    savefig(figh_filtadisp,'FilteredDispForGradient');
                    savefig(figh_exx,'exx_strain');
                    savefig(figh_eyy,'eyy_strain');
                    savefig(figh_exy,'exy_strain');
                    savefig(figh_eyx,'eyx_strain');
                    savefig(figh_gxy,'gxy_pure_shear_strain');
                    savefig(figh_theta,'fixed_body_rotation');
                    savefig(figh_vonMises,'vonMises_strain');
                    savefig(figh_gradcomponentmag,'GradientComponentsMagnitude');
                    savefig(figh_DSCscatter,'StrainFilterDisplacementScatterplot');
                end
                
                % Write a text file giving the conditions 
                mystr = sprintf('This file describes the parameters used when generating strain maps along the following folderpath and foldername:\n');
                mystr = sprintf('%sFolderpath: %s\nFoldername: %s\n\n',mystr,obj.saved_plots_folderpath,obj.saved_plots_foldername);
                if sign_convention(1) == 0 && sign_convention(2) == 0
                    mystr = sprintf('%sMoire rotation has not been removed (rotation plot is total interlayer rotation).\n',mystr);
                else
                    mystr = sprintf('%sMoire rotation has been removed (rotation plot is interlayer rotation due to reconstruction only).\n',mystr);
                    mystr = sprintf('%sSign convention for Moire rotation removal is %d, %d.\n',mystr,sign_convention(1),sign_convention(2));
                end
                if divide_by_two_for_intralayer
                    mystr = sprintf('%sAll maps have been divided by two to correspond to the strain and rotation in an individual layer of graphene.\n',mystr);
                else
                    mystr = sprintf('%sAll maps are left at original measured magnitudes, corresponding to interlayer quantities.\n',mystr);
                end
                mystr = sprintf('%sX-axis has been set perpendicular to the %dth saddle point domains.',mystr,SP_num_for_xaxis);
                mystr = sprintf('%sFilter settings are saved in the corresponding file "filterstruct_used_for_strain_maps.mat".\n',mystr);
                mystr = sprintf('%sPhase-unwrapped displacement maps were trimmed by %d pixels along image boundary before making strain maps.\n',mystr,trim_edge_pixel_num(1));
                if ~isempty(trim_tear_pixel_num)
                    mystr = sprintf('%sPhase-unwrapped displacement maps were trimmed by %d along tear mask before making strain maps.\n',mystr,trim_tear_pixel_num);
                end
                mystr = sprintf('%s\nThis dataset''s filename is: %s\n',mystr,obj.filename);
                mystr = sprintf('%sThis report auto-generated at %s.\n',mystr,datestr(clock));
                
                save('filterstruct_used_for_strain_maps.mat','filterstruct');
                
                fid = fopen('StrainMapParameterSettings.txt','wt');
                fprintf(fid,'%s\n',mystr);
                fclose(fid);
                cd(currentd);
            end
        end
        
        
        
        
        
        
        % This and the corresponding plotting function are the ones that
        % interface into Colin's TGV code. The final and best functions.
        % NPK 06/13/2020
        %
        % pre_filterstruct will allow light filtering of the unwrapped
        % displacement field before TGV. But we don't want too much! Just
        % enough to get some of the outliers out.
        function computeStrainMaps3(obj,trim_value,pre_filterstruct,TGVsettings_id,rotcalibration_ID,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid)
            if (nargin < 3), pre_filterstruct = []; end;
            if (nargin < 2), trim_value = 0; end;
            MAKE_TEST_PLOTS = true;
            obj.trim_value = sum(trim_value);  % So if trim value is more than one element, the maps will be made from the correct trim factor. 
            obj.trim_tear_value = trim_tear_value;
            unwrapped_disp_field_org = obj.annealed_dfield;
            % Edge trimming adjustment not included in Colin's original
            % function. This will allow auto colormap.
            uX = unwrapped_disp_field_org(:,:,1);
            uY = unwrapped_disp_field_org(:,:,2);
            if ~isempty(pre_filterstruct)
                plotflag = false;
                [ uX,uY ] = filterDisplacement( uX,uY,pre_filterstruct,plotflag,obj );
            end
            if trim_value(1) > 0
                uX = trimArray(uX,trim_value(1));
                uY = trimArray(uY,trim_value(1));
            end
            % Convert displacement field into nm.
            uX = 0.1*uX;  % converts the X displacement to nm
            uY = 0.1*uY;  % converts the y displacement to nm
            
            % Only do the mask nanset after the TGV, because TGV cannot
            % handle nans in the image.
% % %             if ~isempty(obj.tear_mask) && ~isempty(trim_tear_value)
% % %                 uX(obj.tear_mask) = 0;
% % %                 uY(obj.tear_mask) = 0;
% % %                 uXorg = uX;
% % %                 uYorg = uY;
% % %             end
            
            unwrapped_disp_field = cat(3,uX,uY);
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            if trim_value(1) > 0
                xbase((end-trim_value(1)+1):end) = [];
                xbase(1:trim_value(1)) = [];
                ybase((end-trim_value(1)+1):end) = [];
                ybase(1:trim_value(1)) = [];
            end
            obj.xbase_for_strainmaps = xbase;
            obj.ybase_for_strainmaps = ybase;
            
            
            
            
            
%             sStrain.flip_x_disp = true;
            sStrain.flip_x_disp = false;  
            % NPK note: Colin originally had flipping of the
            % xdisp field set to true, where everything is moved to
            % negative. I wasn't doing this, so for comparison run under
            % false.
            sStrain.flip_y_disp = false;
            sStrain.flagNormalizeDisp = false;  % Remove mean plane from displacement maps
%             sStrain.strainRange = [-1 1]*0.10;  % for plotting
            % NPK note: will probably want this to go auto for fixed zero
            % colormapping.
            sStrain.strainRange = 'auto';
%             USE_NPK_SETTINGS = true;
            switch TGVsettings_id
                case 'DS1'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS2testing'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 512;  % number of iterations
                    sStrain.TVG_alpha = 0.2;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.05;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS3'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS4'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS5'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS6'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS8'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS10'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 20;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS15'  % This could arguably be 8 iterations as well.
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS18'  % This could arguably be 8 iterations as well.
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS15quant'  % This could arguably be 8 iterations as well.
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 4;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS18quant'  % This could arguably be 8 iterations as well.
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 4;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS26'  % Colin's original
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 64;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS30S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 12;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS8S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS9S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS10S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS11S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS12S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS15S2'  % Might need to build a quant filter for this one too.
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 32;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS3S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 8;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'DS1S2'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 16;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
                case 'minimal'
                    sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
                    sStrain.TVG_iter = 1;  % number of iterations
                    sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
                    sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
                    sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
                    sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
                    sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
            end

%             if USE_NPK_SETTINGS   % Trying to weaken the TVG filter from what Colin originally had
%                 sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
%                 sStrain.TVG_iter = 0;  % number of iterations
%                 sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
%                 sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
%                 sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
%                 sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
%                 sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
%             else
%                 sStrain.TVG_padding = 20; % This parameter is a built-in padArray in the TVG function.
%                 sStrain.TVG_iter = 64;  % number of iterations
%                 sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
%                 sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
%                 sStrain.TVG_theta1 = 5;  % limiting parameter for TVG
%                 sStrain.TVG_theta2 = 5;  % limiting parameter for TVG
%                 sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising
%             end
            % Stepsize adjustment not included in Colin's original
            % function. 
            sStrain.stepsize = obj.scan_stepsize;
            
            [sStrain] = strainCalc01NPK(unwrapped_disp_field,sStrain,MAKE_TEST_PLOTS);
            
            % Only do the mask nanset after the TGV, because TGV cannot
            % handle nans in the image.
            if ~isempty(obj.tear_mask) && ~isempty(trim_tear_value)
                    tear_boundary_mask = boundarymask(obj.tear_mask);
                    tear_boundary_mask_thick = tear_boundary_mask;
                    for i = 1:trim_tear_value
                        tear_boundary_mask_thick = boundarymask(tear_boundary_mask_thick) | tear_boundary_mask_thick;
                        tear_boundary_mask_thick(obj.tear_mask) = false;
                    end
                uX(tear_boundary_mask_thick) = nan;
                uY(tear_boundary_mask_thick) = nan;
                sStrain.xDispTVG(tear_boundary_mask_thick) = nan;
                sStrain.yDispTVG(tear_boundary_mask_thick) = nan;
                sStrain.strainExx(tear_boundary_mask_thick) = nan;
                sStrain.strainExy(tear_boundary_mask_thick) = nan;
                sStrain.strainEyx(tear_boundary_mask_thick) = nan;
                sStrain.strainEyy(tear_boundary_mask_thick) = nan;
            end
            
            obj.StrainStruct = sStrain;
            % Filter again here in case needed, to stop TGV cusping
            if numel(trim_value) > 1
                if trim_value(2) > 0
                    xbase((end-trim_value(2)+1):end) = [];
                    xbase(1:trim_value(2)) = [];
                    ybase((end-trim_value(2)+1):end) = [];
                    ybase(1:trim_value(2)) = [];
                end
                obj.xbase_for_strainmaps = xbase;
                obj.ybase_for_strainmaps = ybase;
                sStrain.xDispTVG = trimArray(sStrain.xDispTVG,trim_value(2));
                sStrain.yDispTVG = trimArray(sStrain.yDispTVG,trim_value(2));
                sStrain.strainExx = trimArray(sStrain.strainExx,trim_value(2));
                sStrain.strainExy = trimArray(sStrain.strainExy,trim_value(2));
                sStrain.strainEyx = trimArray(sStrain.strainEyx,trim_value(2));
                sStrain.strainEyy = trimArray(sStrain.strainEyy,trim_value(2));
                uX = trimArray(uX,trim_value(2));
                uY = trimArray(uY,trim_value(2));
            end
            
            obj.dfield_filtered_for_strain = cat(3,sStrain.xDispTVG,sStrain.yDispTVG);
            % Note that the exx, exy, etc. here are just the raw derivative
            % results (dUinterX_dorgY, e.g.) in the original coordinate
            % system. So we should not update instance variables yet, which
            % assume that the various properties have been satisfied
            % (divide by 2, rotate, whatever is desired).
            
            % Perform rotation of the strain field. This is the same way as
            % under makeStrainMapsFromBlinking2()
            %
            % New 06/30/2020: ability to do the first-principals rotational
            % calibration. rotcalibration_ID is the new variable.
            % Otherwise, calibration needs to be wrt the x-axis (purple
            % soliton wall), and tensor rotation from there on.
            
            switch rotcalibration_ID
                case 'SP'
                    [transformed_coords,rotated_basis,xh1,xh2,yh1,yh2,~,basisrotangle] = obj.changeDisplacementBasis2(1);
                    basisrotangle = rad2deg(basisrotangle);
                case 'nanoparticle'
                    obj.getDiffractionPatternRotation('integration disks');
                    [realspacerotation_deg,SProtation_deg] = obj.getRealSpaceRotationAngle();
                    basisrotangle = realspacerotation_deg;
                    rsrot_rad = deg2rad(basisrotangle);
                    rotated_basis = [cos(rsrot_rad), -sin(rsrot_rad); sin(rsrot_rad), cos(rsrot_rad)];
                case 'nanoparticle mod'
                    obj.getDiffractionPatternRotation('integration disks');
                    [realspacerotation_deg,SProtation_deg] = obj.getRealSpaceRotationAngle();
                    basisrotangle = mod(realspacerotation_deg,180);
                    rsrot_rad = deg2rad(basisrotangle);
                    rotated_basis = [cos(rsrot_rad), -sin(rsrot_rad); sin(rsrot_rad), cos(rsrot_rad)];
                otherwise
                        error('Please choose a valid option for rotational calibration.');
            end
            
            obj.rotated_basis = rotated_basis;
            obj.initial_tensor_angle = basisrotangle;
            % Retrieve the gradient of the filtered displacement data from sStrain.
            dUinterX_dorgx = sStrain.strainExx;
            dUinterX_dorgy = sStrain.strainExy;
            dUinterY_dorgx = sStrain.strainEyx;
            dUinterY_dorgy = sStrain.strainEyy;
            % Next, compute the directional derivatives of the x and y
            % components of the interlayer strain. Do this by dot
            % products with a unit normal vector, using the directional
            % derivative formula.
            xdir = rotated_basis(:,1)';
            ydir = rotated_basis(:,2)';
            dUinterX = vertcat(dUinterX_dorgx(:)',dUinterX_dorgy(:)');
            dUinterX_dx = reshape((xdir*dUinterX)',size(uX));
            dUinterX_dy = reshape((ydir*dUinterX)',size(uX));
            dUinterY = vertcat(dUinterY_dorgx(:)',dUinterY_dorgy(:)');
            dUinterY_dx = reshape((xdir*dUinterY)',size(uX));
            dUinterY_dy = reshape((ydir*dUinterY)',size(uX));
            % Now, these derivatives represent the interlayer strain in
            % the new coordinate system. But we want to convert them to
            % intralayer strain using the moire twist angle, estimated
            % independently.
            if ~isempty(obj.moire_angle_estimates)
                mangles = trimArray(obj.moire_angle_estimates,trim_value(1));  % Will need to crop this to match.
                if numel(trim_value) > 1
                    mangles = trimArray(obj.moire_angle_estimates,trim_value(2));  % Will need to crop this to match.
                end
            else
                warning('There is no Moire angle set. This could lead to inaccurate simple shear strain maps.');
                tf = input('Do you wish to proceed? 1/0:');
                if ~tf
                    error('Aborting program');
                end
                mangles = zeros(size(dUinterX_dx));
            end
            
            if divide_by_two_for_intralayer
                obj.exx = 0.5*(dUinterX_dx);
                obj.exy = 0.5*(dUinterX_dy + sign_convention(1)*deg2rad(mangles));
                obj.eyx = 0.5*(dUinterY_dx + sign_convention(2)*deg2rad(mangles));
                obj.eyy = 0.5*(dUinterY_dy);
            else
                obj.exx = (dUinterX_dx);
                obj.exy = (dUinterX_dy + sign_convention(1)*deg2rad(mangles));
                obj.eyx = (dUinterY_dx + sign_convention(2)*deg2rad(mangles));
                obj.eyy = (dUinterY_dy);
            end
            % First, calculate pure shear and FBR, because we will use them
            % when performing the tensor rotations.
            %
            % 07/16/2020: Notation alert! What we originally called
            % "obj.gxy" here is not actually the engineering pure shear
            % strain, but rather the tensorial pure shear strain. 
            obj.gxy = 0.5*(obj.exy + obj.eyx);
            obj.fixed_body_rotation = 0.5*(obj.eyx - obj.exy);
            obj.FBR_type = 'radians';
            % 05/30/2020: reworked handling of the cross derivative
            % terms after figuring out the issue with the vonMises
            % strain invariant.
            %
            % Note that these values are just reparameterizations of
            % the simple shear strain. They are not to be multiplied by
            % 1/2 again; that is simply part of the definition for the
            % strain tensor.
            
            % Here we can perform tensor rotations to reference the strain
            % against other soliton walls. This will, say, help us figure
            % out the orientation of the principal strain through the soliton walls
            % for distorted samples. It will presumably not be the same for
            % all soliton walls if the lattce is non-hexagonal.
            %
            % However, first we still have to do the directional derivative
            % orientation to produce a strain tensor.
            %             [sStrain] = strainCalc02NPK(sStrain,-basisrotangle);
            if ~isempty(tensor_rotate_SPid) % Do the tensor rotation
                old_tensor_angle = deg2rad(basisrotangle);
                % note: this function demands that 
                [~,~,~,~,~,~,~,basisrotangle2] = obj.changeDisplacementBasis2(tensor_rotate_SPid);
                new_tensor_angle = basisrotangle2;
                [exx_rot,exy_rot,eyx_rot,eyy_rot,gxy_rot,newbasis] = obj.rotateStrainTensor(old_tensor_angle,new_tensor_angle,obj.rotated_basis);
                obj.rotated_basis = newbasis;
                obj.initial_tensor_angle = rad2deg(new_tensor_angle);
                obj.exx = exx_rot;
                obj.exy = exy_rot;
                obj.eyx = eyx_rot;
                obj.eyy = eyy_rot;
                obj.gxy = gxy_rot;
            end
            
            
            obj.exx = 100*obj.exx;
            obj.exy = 100*obj.exy;
            obj.eyx = 100*obj.eyx;
            obj.eyy = 100*obj.eyy;
            obj.gxy = 100*obj.gxy;
            % 05/30/2020: express fixed-body rotation in terms of degrees.
            obj.fixed_body_rotation = rad2deg(obj.fixed_body_rotation);
            obj.FBR_type = 'degrees';
            obj.vonMises = obj.gxy.^2 - obj.exx.*obj.eyy;
            
            % Use Colin's third function to obtain the principal strains.
            [imageRGB,E1,E2,theta,Ediff,Edil] = strainCalc03NPK(obj.exx,obj.eyy,obj.gxy);
            % Store all of these; may be useful.
            obj.principal_strain_1 = E1;  % Emax
            obj.principal_strain_2 = E2;  % Emin
            obj.principal_angle = rad2deg(theta);
%             obj.principal_strain_max = Emax;
            obj.principal_strain_diff = Ediff;
            obj.principal_strain_dilatation = Edil;
            obj.principal_strain_3Dcolor = imageRGB;
            
            obj.sign_convention = sign_convention;  % store for plots
            obj.divide_by_two_for_intralayer = divide_by_two_for_intralayer;
            
            
            
            
            % Plug in the results of Colin's functions to the instance
            % variable slots.
            %             obj.exx = sStrain.strainExx;
            %             obj.exy = sStrain.strainExy;
            %             obj.eyx = sStrain.strainEyx;
            %             obj.eyy = sStrain.strainEyy;
        end
        
        
        
        
        
        
        
        
        % This function requires that computeStrainMaps3 have already been
        % run. All relevant information is stored.
        %
        % Modified 07/01/2020 to try to get the data tips to work with the
        % overlays.
        function plotStrainMaps3(obj,saveplotflag,overlay_registration,colormaptype,crop_range)
            trim_edge_pixel_num = obj.trim_value;
            sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
            if nargin < 5
                crop_range = [];
            end
            
            if ~isempty(crop_range)
                obj.exx = obj.exx(crop_range{1},crop_range{2});
                obj.eyy = obj.eyy(crop_range{1},crop_range{2});
                obj.exy = obj.exy(crop_range{1},crop_range{2});
                obj.eyx = obj.eyx(crop_range{1},crop_range{2});
                obj.gxy = obj.gxy(crop_range{1},crop_range{2});
                obj.fixed_body_rotation = obj.fixed_body_rotation(crop_range{1},crop_range{2});
                obj.principal_strain_1 = obj.principal_strain_1(crop_range{1},crop_range{2});
                obj.principal_strain_2 = obj.principal_strain_2(crop_range{1},crop_range{2});
                obj.dfield_filtered_for_strain = obj.dfield_filtered_for_strain(crop_range{1},crop_range{2},:);
%                 obj.principal_strain_max = obj.principal_strain_max(crop_range{1},crop_range{2});
                obj.principal_strain_diff = obj.principal_strain_diff(crop_range{1},crop_range{2});
                obj.principal_strain_dilatation = obj.principal_strain_dilatation(crop_range{1},crop_range{2});
                obj.xbase_for_strainmaps = obj.xbase_for_strainmaps(crop_range{2});
                obj.ybase_for_strainmaps = obj.ybase_for_strainmaps(crop_range{1});
%                 eyy
%                 exy
%                 eyx
%                 gxy
%                 fixed_body_rotation
%                 dfield_filtered_for_strain
%                 principal_strain_1
%                 principal_strain_2
%                 principal_strain_max
%                 principal_strain_diff
%                 principal_strain_dilatation
                sw = sw(crop_range{1},crop_range{2});
            end
            
            ALPHA_CONSTANT = 0.25;
            ALPHA_CONSTANT2 = 0; % Simply trying to trick the data tips into working.
            if isa(colormaptype,'cell')
                colormaptype_principal = colormaptype{2};
                colormaptype = colormaptype{1};
                scale_principal_color_map_flag = false;
            else
                colormaptype_principal = colormaptype;
                scale_principal_color_map_flag = true;
            end
            
            divide_by_two_for_intralayer = obj.divide_by_two_for_intralayer;
            sign_convention = obj.sign_convention;
            rotated_basis = obj.rotated_basis;
            
            shadingtype = 'flat';
            scale_color_map_flag = true;
            
            xbase = obj.xbase_for_strainmaps;
            ybase = obj.ybase_for_strainmaps;
            uX = obj.dfield_filtered_for_strain(:,:,1);
            uY = obj.dfield_filtered_for_strain(:,:,2);
            % Factor of ten needed here because of the reduced zone
            % conversion
            [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( 10*[uX(:),uY(:)] );
            uXr = reshape(reduced_zone_disps(:,1),size(uX));
            uYr = reshape(reduced_zone_disps(:,2),size(uX));
            [figh_filtadisp,~] = makeOutsideCustomDisplacementColorPlot(obj,uXr,uYr);
            title('Reduced zone representation of filtered data for strain mapping');
            
            
            figh_DSCscatter = figure;
            DSC1 = uX;
            DSC2 = uY;
            scatter(DSC1(:)*10,DSC2(:)*10,3,'Filled');
            axh = gca;
            plotFullDisplacementHexagons(axh);
%             xlim([-1.5,1.5]);
%             ylim([-0.1,1.5]);
            xlabel('Cartesian displacement component 1');
            ylabel('Cartesian displacement component 2');
            title('Displacement vector scatterplot after pre-strain filter');
            
            
            xcolor = 'g';
            ycolor = 'm';
            
            
            %%%%% exx
            figh_exx = figure;
            plot_exx = obj.exx; % trimArray(obj.exx,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_exx);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_exx,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('exx intralayer normal strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('exx interlayer normal strain');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
                
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_exx);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            %%%%% eyy
            figh_eyy = figure;
            plot_eyy = obj.eyy;%trimArray(obj.eyy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_eyy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_eyy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('eyy intralayer normal strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('eyy interlayer normal strain');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_eyy);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            %%%%% exy figure
            figh_exy = figure;
            plot_exy = obj.exy; %trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_exy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_exy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('exy intralayer simple shear strain (cross derivative)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('exy interlayer simple shear strain (cross derivative)');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_exy);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            %%%%% eyx figure
            % Note: in the eventual nomenclature of the TBG paper, these
            % are all syx sxy simple shear strains. 
            figh_eyx = figure;
            plot_eyx = obj.eyx;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_eyx);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_eyx,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('eyx intralayer simple shear strain (cross derivative)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('eyx interlayer simple shear strain (cross derivative)');
                ylabel(cbh,'interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_eyx);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            %%%%% vonMises strain
            vM = obj.vonMises;
            figh_vonMises = figure;
            imagesc(xbase,ybase,vM);
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('Intralayer vonMises Strain Invariant (fixed body rotation removed)');
                ylabel(cbh,'(Intralayer strain %)^2');
            else
                title('Interlayer vonMises Strain Invariant (fixed body rotation removed)');
                ylabel(cbh,'(Interlayer strain %)^2');
            end
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,vM);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            
            %%%% gxy tensorial pure shear strain plot
            figh_gxy = figure;
            plot_gxy = obj.gxy;%trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_gxy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_gxy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('gyx intralayer tensorial pure shear strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('gyx interlayer tensorial pure shear strain');
                ylabel(cbh,'Interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_gxy);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            %%%% gammaxy engineering pure shear strain plot
            figh_gammaxy = figure;
            plot_gxy = 2*obj.gxy;%trimArray(obj.gxy,trim_edge_pixel_num);
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_gxy);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_gxy,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype)
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('gammayx intralayer engineering pure shear strain');
                ylabel(cbh,'Intralayer strain %');
            else
                title('gammayx interlayer engineering pure shear strain');
                ylabel(cbh,'Interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_gxy);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            %%%% Fixed body rotation figure
            figh_theta = figure;
            plot_theta = obj.fixed_body_rotation;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_theta);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_theta,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if sign_convention(1) == 0 && sign_convention(2) == 0
                if divide_by_two_for_intralayer
                    title('Intralayer fixed body rotation, Moire rotation included');
                    ylabel(cbh,'Intralayer total rotation (degrees)');
                else
                    title('Interlayer fixed body rotation, Moire rotation included');
                    ylabel(cbh,'Interlayer total rotation (degrees)');
                end
            else
                if divide_by_two_for_intralayer
                    title('Intralayer fixed body rotation, Moire rotation removed');
                    ylabel(cbh,'Intralayer reconstruction rotation (degrees)');
                else
                    title('Interlayer fixed body rotation, Moire rotation removed');
                    ylabel(cbh,'Interlayer reconstruction rotation (degrees)');
                end
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_theta);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            
            %%% Principal Strain 1 (max)
            figh_E1 = figure;
            plot_E1 = obj.principal_strain_1;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_E1);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_E1,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_principal_color_map_flag
                scaleColorMap( colormaptype_principal, 0);
            else
                colormap(colormaptype_principal);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('Principal strain max component (intralayer)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('Principal strain max component (interlayer)');
                ylabel(cbh,'Interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_E1);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            
            %%% Principal Strain 2 (min)
            figh_E2 = figure;
            plot_E2 = obj.principal_strain_2;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_E2);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_E2,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_principal_color_map_flag
                scaleColorMap( colormaptype_principal, 0);
            else
                colormap(colormaptype_principal);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('Principal strain min component (intralayer)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('Principal strain min component (interlayer)');
                ylabel(cbh,'Interlayer strain %');
            end
            % set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT+0.3;
                h2 = imagesc(xbase,ybase,plot_E2);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            
            %%% Principal Angle
            figh_Eangle = figure;
            plot_Eangle = obj.principal_angle;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_Eangle);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_Eangle,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            colormap(hsv);
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            title('Principal angle');
            ylabel(cbh,'Degrees');
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,cat(3,sw,sw,sw));
                h.AlphaData = 0.25;
                h2 = imagesc(xbase,ybase,plot_Eangle);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            
            % This plot unneeded after the angle unwrapping, making it
            % explicit that E1 and E2 are Emax and Emin.
% % %             %%% Principal Strain Max
% % %             figh_Emax = figure;
% % %             plot_Emax = obj.principal_strain_max;
% % %             if strcmp(shadingtype,'flat')
% % %                 imagesc(xbase,ybase,plot_Emax);
% % %             elseif strcmp(shadingtype,'interp')
% % %                 contourf(xbase,ybase,plot_Emax,NUM_CONTOUR_LINES,'LineStyle','none');
% % %                 shading(gca,shadingtype);
% % %             end
% % %             if scale_principal_color_map_flag
% % %                 scaleColorMap( colormaptype_principal, 0);
% % %             else
% % %                 colormap(colormaptype_principal);
% % %             end
% % %             axis equal
% % %             xlim([xbase(1) xbase(end)]);
% % %             ylim([ybase(1) ybase(end)]);
% % %             cbh = colorbar;
% % %             if divide_by_two_for_intralayer
% % %                 title('Principal strain max component (intralayer)');
% % %                 ylabel(cbh,'Intralayer strain %');
% % %             else
% % %                 title('Principal strain max component (interlayer)');
% % %                 ylabel(cbh,'Interlayer strain %');
% % %             end
% % %             xlabel('x (nm)');
% % %             ylabel('y (nm)');
% % %             set(gca,'ydir','normal');
% % %             hold on
% % %             xl = get(gca,'xlim');
% % %             yl = get(gca,'ylim');
% % %             qcoords = [xl(1),yl(1)] + 10;
% % %             quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
% % %             quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
% % %             if overlay_registration
% % %                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
% % %                 h = imagesc(xbase,ybase,sw);
% % %                 h.AlphaData = ALPHA_CONSTANT;
% % %             end
% % %             
            
            
            %%% Principal Strain Quiver
            % NOTE: check if should be Emax or E1.
            % Resolved: they are the same
            figh_principalquiver = figure;
            try
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                hold on
            end
            % These are stored in degrees, but need rad for the unit vector
            % calculation.
            %
            % Need to add initial tensor angle because tensor rotation is
            % referenced on top of the original coordinate orientation of
            % the tensor.
            angles = deg2rad(obj.principal_angle(:)) + obj.initial_tensor_angle;
%             angles = deg2rad(obj.principal_angle(:));
            unit_vecs = horzcat(cos(angles),sin(angles));
            vecs = unit_vecs.*obj.principal_strain_1(:);
            [xspace,yspace] = meshgrid(xbase,ybase);
            qh1 = quiver(xspace(:),yspace(:),vecs(:,1),vecs(:,2),'k');
            hold on;
            qh2 = quiver(xspace(:),yspace(:),-vecs(:,1),-vecs(:,2),'k');
            qh1.ShowArrowHead = 'off';
            qh2.ShowArrowHead = 'off';
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            axis equal
            title('Principal strain quiverplot (full)');
            end
            
            
            figh_reducedprincipalquiver = figure;
            try
            scaling = 0.4;
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                hold on
            end
            xspace_full = reshape(xspace,size(obj.exx));
            yspace_full = reshape(yspace,size(obj.exx));
            veccomp1 = reshape(vecs(:,1),size(obj.exx));
            veccomp2 = reshape(vecs(:,2),size(obj.exx));
            [indend1,indend2] = size(obj.exx);
            downsample_factor = 6;
            ind1rng = 1:downsample_factor:indend1;
            ind2rng = 1:downsample_factor:indend2;
            xspace_red = xspace_full(ind1rng,ind2rng);
            yspace_red = yspace_full(ind1rng,ind2rng);
            veccomp1red = veccomp1(ind1rng,ind2rng);
            veccomp2red = veccomp2(ind1rng,ind2rng);
            qh1 = quiver(xspace_red(:),yspace_red(:),veccomp1red(:),veccomp2red(:),scaling,'k');
            hold on;
            qh2 = quiver(xspace_red(:),yspace_red(:),-veccomp1red(:),-veccomp2red(:),scaling,'k');
            qh1.ShowArrowHead = 'off';
            qh2.ShowArrowHead = 'off';
            qh1.LineWidth = 1.0;
            qh2.LineWidth = 1.0;
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            axis equal
            title('Principal strain quiverplot (downsampled)');
            end
            
            
            %%% Principal Strain Difference
            figh_Ediff = figure;
            plot_Ediff = obj.principal_strain_diff;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_Ediff);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_Ediff,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_principal_color_map_flag
                scaleColorMap( colormaptype_principal, 0);
            else
                colormap(colormaptype_principal);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('Principal strain component difference (intralayer)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('Principal strain component difference (interlayer)');
                ylabel(cbh,'Interlayer strain %');
            end
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_Ediff);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            
            
            %%% Principal Strain Dilatation
            % This one is always informative in the same colormap as the
            % regular strain plots.
            figh_Edil = figure;
            plot_Edil = obj.principal_strain_dilatation;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_Edil);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_Edil,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            if scale_color_map_flag
                scaleColorMap( colormaptype, 0);
            else
                colormap(colormaptype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            cbh = colorbar;
            if divide_by_two_for_intralayer
                title('Dilatation strain (intralayer)');
                ylabel(cbh,'Intralayer strain %');
            else
                title('Dilatation strain (interlayer)');
                ylabel(cbh,'Interlayer strain %');
            end
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
                h2 = imagesc(xbase,ybase,plot_Edil);
                h2.AlphaData = ALPHA_CONSTANT2;
            end
            
            
            %%% Principal Strain 3D colorization
            figh_Ecolor = figure;
            plot_Ecolor = obj.principal_strain_3Dcolor;
            if strcmp(shadingtype,'flat')
                imagesc(xbase,ybase,plot_Ecolor);
            elseif strcmp(shadingtype,'interp')
                contourf(xbase,ybase,plot_Ecolor,NUM_CONTOUR_LINES,'LineStyle','none');
                shading(gca,shadingtype);
            end
            axis equal
            xlim([xbase(1) xbase(end)]);
            ylim([ybase(1) ybase(end)]);
            if divide_by_two_for_intralayer
                title('3D colorized principal strains (intralayer)');
            else
                title('3D colorized principal strains (interlayer)');
            end
            xlabel('x (nm)');
            ylabel('y (nm)');
            set(gca,'ydir','normal');
            hold on
            xl = get(gca,'xlim');
            yl = get(gca,'ylim');
            qcoords = [xl(1),yl(1)] + 10;
            quiver(qcoords(1),qcoords(2),rotated_basis(1,1)*10,rotated_basis(2,1)*10,'Color',xcolor);
            quiver(qcoords(1),qcoords(2),rotated_basis(1,2)*10,rotated_basis(2,2)*10,'Color',ycolor);
            if overlay_registration
%                 sw = trimArray(full(obj.soliton_walls_merged),trim_edge_pixel_num(1));
                h = imagesc(xbase,ybase,sw);
                h.AlphaData = ALPHA_CONSTANT;
            end
        
        
            
            %%%% Make statistics
            eyx_percentiles = prctile(obj.eyx(:),0:100);
            gradientsum_percentiles = prctile(obj.NPK_strain_mag(:),0:100);
            eyx_values = obj.eyx;
            gradientsum_values = obj.NPK_strain_mag;
            
            if saveplotflag
                currentd = pwd;
                cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                
                if overlay_registration
                    mkdir('OverlayStrainFigures');
                    cd('OverlayStrainFigures');
                    saveas(figh_filtadisp,'FilteredDispForGradient_overlay.png');
                    saveas(figh_exx,'exx_strain_overlay.png');
                    saveas(figh_eyy,'eyy_strain_overlay.png');
                    saveas(figh_exy,'exy_strain_overlay.png');
                    saveas(figh_eyx,'eyx_strain_overlay.png');
                    saveas(figh_gxy,'gxy_tensorial_pure_shear_strain_overlay.png');
                    saveas(figh_gammaxy,'gxy_engineering_pure_shear_strain_overlay.png');
                    saveas(figh_theta,'fixed_body_rotation_overlay.png');
                    saveas(figh_vonMises,'vonMises_strain_overlay.png');
%                     saveas(figh_gradcomponentmag,'GradientComponentsMagnitude_overlay.png');
                    saveas(figh_DSCscatter,'StrainFilterDisplacementScatterplot.png');
                    saveas(figh_Ecolor,'PrincipalStrain3DColorization_overlay.png');
                    saveas(figh_Edil,'PrincipalStrainDilatation_overlay.png');
                    saveas(figh_Ediff,'PrincipalStrainE1E2Difference_overlay.png');
                    saveas(figh_reducedprincipalquiver,'PrincipalStrainReducedQuiver_overlay.png');
                    saveas(figh_principalquiver,'PrincipalStrainFullQuiver_overlay.png');
                    saveas(figh_Eangle,'PrincipalAngle_overlay.png');
                    saveas(figh_E1,'PrincipalStrainMaxComponent_overlay.png');
                    saveas(figh_E2,'PrincipalStrainMinComponent_overlay.png');                    
                    
                    savefig(figh_filtadisp,'FilteredDispForGradient_overlay');
                    savefig(figh_exx,'exx_strain_overlay');
                    savefig(figh_eyy,'eyy_strain_overlay');
                    savefig(figh_exy,'exy_strain_overlay');
                    savefig(figh_eyx,'eyx_strain_overlay');
                    savefig(figh_gxy,'gxy_tensorial_pure_shear_strain_overlay');
                    savefig(figh_gammaxy,'gxy_engineering_pure_shear_strain_overlay');
                    savefig(figh_theta,'fixed_body_rotation_overlay');
                    savefig(figh_vonMises,'vonMises_strain_overlay');
%                     savefig(figh_gradcomponentmag,'GradientComponentsMagnitude_overlay');
                    savefig(figh_DSCscatter,'StrainFilterDisplacementScatterplot');
                    savefig(figh_Ecolor,'PrincipalStrain3DColorization_overlay');
                    savefig(figh_Edil,'PrincipalStrainDilatation_overlay');
                    savefig(figh_Ediff,'PrincipalStrainE1E2Difference_overlay');
                    savefig(figh_reducedprincipalquiver,'PrincipalStrainReducedQuiver_overlay');
                    savefig(figh_principalquiver,'PrincipalStrainFullQuiver_overlay');
                    savefig(figh_Eangle,'PrincipalAngle_overlay');
                    savefig(figh_E1,'PrincipalStrainMaxComponent_overlay');
                    savefig(figh_E2,'PrincipalStrainMinComponent_overlay');  
                else
                    mkdir('StrainFigures');
                    cd('StrainFigures');
                    saveas(figh_filtadisp,'FilteredDispForGradient.png');
                    saveas(figh_exx,'exx_strain.png');
                    saveas(figh_eyy,'eyy_strain.png');
                    saveas(figh_exy,'exy_strain.png');
                    saveas(figh_eyx,'eyx_strain.png');
                    saveas(figh_gxy,'gxy_tensorial_pure_shear_strain.png');
                    saveas(figh_gammaxy,'gxy_engineering_pure_shear_strain.png');
                    saveas(figh_theta,'fixed_body_rotation.png');
                    saveas(figh_vonMises,'vonMises_strain.png');
%                     saveas(figh_gradcomponentmag,'GradientComponentsMagnitude.png');
                    saveas(figh_DSCscatter,'StrainFilterDisplacementScatterplot.png');
                    saveas(figh_Ecolor,'PrincipalStrain3DColorization.png');
                    saveas(figh_Edil,'PrincipalStrainDilatation.png');
                    saveas(figh_Ediff,'PrincipalStrainE1E2Difference.png');
                    saveas(figh_reducedprincipalquiver,'PrincipalStrainReducedQuiver.png');
                    saveas(figh_principalquiver,'PrincipalStrainFullQuiver.png');
                    saveas(figh_Eangle,'PrincipalAngle.png');
                    saveas(figh_E1,'PrincipalStrainMaxComponent.png');
                    saveas(figh_E2,'PrincipalStrainMinComponent.png'); 
                    
                    savefig(figh_filtadisp,'FilteredDispForGradient');
                    savefig(figh_exx,'exx_strain');
                    savefig(figh_eyy,'eyy_strain');
                    savefig(figh_exy,'exy_strain');
                    savefig(figh_eyx,'eyx_strain');
                    savefig(figh_gxy,'gxy_tensorial_pure_shear_strain');
                    savefig(figh_gammaxy,'gxy_engineering_pure_shear_strain');
                    savefig(figh_theta,'fixed_body_rotation');
                    savefig(figh_vonMises,'vonMises_strain');
%                     savefig(figh_gradcomponentmag,'GradientComponentsMagnitude');
                    savefig(figh_DSCscatter,'StrainFilterDisplacementScatterplot');
                    savefig(figh_Ecolor,'PrincipalStrain3DColorization');
                    savefig(figh_Edil,'PrincipalStrainDilatation');
                    savefig(figh_Ediff,'PrincipalStrainE1E2Difference');
                    savefig(figh_reducedprincipalquiver,'PrincipalStrainReducedQuiver');
                    savefig(figh_principalquiver,'PrincipalStrainFullQuiver');
                    savefig(figh_Eangle,'PrincipalAngle_overlay');
                    savefig(figh_E1,'PrincipalStrainMaxComponent');
                    savefig(figh_E2,'PrincipalStrainMinComponent'); 
                end
                
                % Write a text file giving the conditions 
                mystr = sprintf('This file describes the parameters used when generating strain maps along the following folderpath and foldername:\n');
                mystr = sprintf('%sFolderpath: %s\nFoldername: %s\n\n',mystr,obj.saved_plots_folderpath,obj.saved_plots_foldername);
                if sign_convention(1) == 0 && sign_convention(2) == 0
                    mystr = sprintf('%sMoire rotation has not been removed (rotation plot is total interlayer rotation).\n',mystr);
                else
                    mystr = sprintf('%sMoire rotation has been removed (rotation plot is interlayer rotation due to reconstruction only).\n',mystr);
                    mystr = sprintf('%sSign convention for Moire rotation removal is %d, %d.\n',mystr,sign_convention(1),sign_convention(2));
                end
                if divide_by_two_for_intralayer
                    mystr = sprintf('%sAll maps have been divided by two to correspond to the strain and rotation in an individual layer of graphene.\n',mystr);
                else
                    mystr = sprintf('%sAll maps are left at original measured magnitudes, corresponding to interlayer quantities.\n',mystr);
                end
%                 mystr = sprintf('%sX-axis has been set perpendicular to the %dth saddle point domains.',mystr,SP_num_for_xaxis);
%                 mystr = sprintf('%sFilter settings are saved in the corresponding file "filterstruct_used_for_strain_maps.mat".\n',mystr);
                mystr = sprintf('%sPhase-unwrapped displacement maps were trimmed by %d pixels along image boundary before making strain maps.\n',mystr,trim_edge_pixel_num(1));
                if ~isempty(obj.trim_tear_value)
                    mystr = sprintf('%sPhase-unwrapped displacement maps were trimmed by %d along tear mask before making strain maps.\n',mystr,obj.trim_tear_value);
                end
                mystr = sprintf('%s\nThis dataset''s filename is: %s\n',mystr,obj.filename);
                mystr = sprintf('%sThis report auto-generated at %s.\n',mystr,datestr(clock));
                
%                 save('filterstruct_used_for_strain_maps.mat','filterstruct');
                
                fid = fopen('StrainMapParameterSettings.txt','wt');
                fprintf(fid,'%s\n',mystr);
                fclose(fid);
                cd(currentd);
            end
        end
        
        
        
        
        
        
        function calculatePrincipalSimpleShearPlot(obj)
            % We want the Moire rotation removed here, because the simple
            % shear argument is saying something about the intralayer
            % strain mechanics.
            warning('This function will return different results depending on whether the Moire rotation has been removed or not.');
            n_pixels = numel(obj.exx);  % has already been trimmed
            obj.shear_axes = zeros(n_pixels,1);
            init_locs = linspace(0.01,pi-0.01,4);
            lb = 0;
            ub = pi;
            options = optimoptions('lsqnonlin');
            options.Display = 'off';
            options.FunctionTolerance = 1e-2;
            options.StepTolerance = 1e-2;
            for i = 1:n_pixels
                if mod(i,100) == 0
                    fprintf('Fitting pixel %d of %d.\n',i,n_pixels);
                end
                this_tensor = [obj.exx(i),obj.exy(i);obj.eyx(i),obj.eyy(i)];
                objfcn = @(theta) simpleShearObjFcn( theta, this_tensor );
                residual_stor = zeros(numel(init_locs),1);
                theta_stor = zeros(numel(init_locs),1);
                for q = 1:numel(init_locs)
                    this_init_theta = init_locs(q);
                    calc_theta = lsqnonlin(objfcn,this_init_theta,lb,ub,options);
                    residual_normal = simpleShearObjFcn( calc_theta, this_tensor );
                    theta_stor(q) = calc_theta;
                    residual_stor(q) = sqrt(sum(residual_normal.^2));
                end
                [~,idx] = min(residual_stor);
                obj.shear_axes(i) = theta_stor(idx);
            end
            obj.shear_axes = reshape(obj.shear_axes,size(obj.exx));
        end
        
        
        
        
        
        
        function makePrincipalSimpleShearPlot(obj)
            % Regenerate the full strain tensor at each location
            exx_rotated = zeros(size(obj.exx));
            exy_rotated = zeros(size(obj.exx));
            eyx_rotated = zeros(size(obj.exx));
            eyy_rotated = zeros(size(obj.exx));
            shear_angles = zeros(size(obj.exx));
            n = numel(obj.shear_axes);
            for i = 1:n
                if mod(i,1000) == 0
                    fprintf('Obtaining rotated tensor for pixel %d of %d.\n',i,n);
                end
                this_tensor = [obj.exx(i),obj.exy(i);obj.eyx(i),obj.eyy(i)];
                this_theta = obj.shear_axes(i);
                rottensor = simpleShearPredfun(this_theta,this_tensor);
                rottensor90 = simpleShearPredfun(this_theta+pi/2,this_tensor);
                % Empirical testing shows that again we only need rotation
                % by 90 degrees, not 360.
                % convention: always put this as larger in (2,1)
%                 if abs(rottensor(2,1)) >= abs(rottensor(1,2))
                if rottensor(2,1) >= rottensor90(2,1)
                    shear_angles(i) = this_theta; 
                    exx_rotated(i) = rottensor(1,1);
                    exy_rotated(i) = rottensor(1,2);
                    eyx_rotated(i) = rottensor(2,1);
                    eyy_rotated(i) = rottensor(2,2);
%                 elseif abs(rottensor90(2,1)) >= abs(rottensor90(1,2))
                elseif true
                    shear_angles(i) = this_theta + pi/2;
                    exx_rotated(i) = rottensor90(1,1);
                    exy_rotated(i) = rottensor90(1,2);
                    eyx_rotated(i) = rottensor90(2,1);
                    eyy_rotated(i) = rottensor90(2,2);
                else
                    error('Violated assumption');
                end
                shear_angles = mod(shear_angles,pi);
            end
            
            figure; 
            imagesc(shear_angles); set(gca,'yDir','normal'); axis equal; colormap(hsv); colorbar;
            title('Shear angle');
            figure; 
            imagesc(exx_rotated); set(gca,'yDir','normal'); axis equal; colormap(fire); colorbar;
            title('exx rotated');
            figure; 
            imagesc(exy_rotated); set(gca,'yDir','normal'); axis equal; colormap(fire); colorbar;
            title('exy rotated');
            figure; 
            imagesc(eyx_rotated); set(gca,'yDir','normal'); axis equal; colormap(fire); colorbar;
            title('eyx rotated');
            figure; 
            imagesc(eyy_rotated); set(gca,'yDir','normal'); axis equal; colormap(fire); colorbar;
            title('eyy rotated');
            
            % Make a quiver plot analogous to principal strain quiver plot
            % (1) Entirely on the eyx component, showing the direction of
            % the pulling. This will require rotating the angle by 90
            % degrees, because it is the y-axis we care about (direction of
            % pulling).
            pull_angle_1 = shear_angles + obj.initial_tensor_angle + pi/2;
            unit_vectors = horzcat(cos(pull_angle_1(:)),sin(pull_angle_1(:)));
            vecs = unit_vectors.*eyx_rotated(:);  % not worrying about sign or direction here, yet.
            figure;
            xbase = obj.xbase_for_strainmaps;
            ybase = obj.ybase_for_strainmaps;
            [xspace,yspace] = meshgrid(xbase,ybase);
            scaling = 0.4;
            qh1 = quiver(xspace(:),yspace(:),vecs(:,1),vecs(:,2),scaling,'k');
            hold on
            qh2 = quiver(xspace(:),yspace(:),-vecs(:,1),-vecs(:,2),scaling,'k');
%             qh1.ShowArrowHead = 'off';
%             qh2.ShowArrowHead = 'off';
            set(qh1,'ShowArrowHead','off');
            set(qh2,'ShowArrowHead','off');
            set(qh1,'LineWidth',1.5);
            set(qh2,'LineWidth',1.5);
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Principal simple shear quiver 1: rotated eyx only');
            axis equal
            
            % (2) Also including the exy component, which has a degree of
            % shearing that is along the x axis
            figure
            pull_angle_2 = shear_angles + obj.initial_tensor_angle;
            unit_vectors = horzcat(cos(pull_angle_2(:)),sin(pull_angle_2(:)));
            newvecs = unit_vectors.*exy_rotated(:);  % not worrying about sign or direction here, yet.
            qh1 = quiver(xspace(:),yspace(:),newvecs(:,1),newvecs(:,2),scaling,'k');
            hold on
            qh2 = quiver(xspace(:),yspace(:),-newvecs(:,1),-newvecs(:,2),scaling,'k');
%             qh1.ShowArrowHead = 'off';
%             qh2.ShowArrowHead = 'off';
            set(qh1,'ShowArrowHead','off');
            set(qh2,'ShowArrowHead','off');
            set(qh1,'LineWidth',1.5);
            set(qh2,'LineWidth',1.5);
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Principal simple shear quiver 2: rotated exy only');
            axis equal
            
            % (3) Combination of the two vectors
            sumvecs = newvecs + vecs;
            figure;
            qh1 = quiver(xspace(:),yspace(:),sumvecs(:,1),sumvecs(:,2),scaling,'k');
            hold on
            qh2 = quiver(xspace(:),yspace(:),-sumvecs(:,1),-sumvecs(:,2),scaling,'k');
%             qh1.ShowArrowHead = 'off';
%             qh2.ShowArrowHead = 'off';
            set(qh1,'ShowArrowHead','off');
            set(qh2,'ShowArrowHead','off');
            set(qh1,'LineWidth',1.5);
            set(qh2,'LineWidth',1.5);
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Principal simple shear quiver 3: Summed simple shear axes');
            axis equal
            
            % (4) Plot showing vectors of both orientations
            figure;
            AB_shear_vectors = vecs;
            AA_shear_vectors = newvecs;
            qh1 = quiver(xspace(:),yspace(:),AB_shear_vectors(:,1),AB_shear_vectors(:,2),scaling,'r');
            hold on
%             qh2 = quiver(xspace(:),yspace(:),-AB_shear_vectors(:,1),-AB_shear_vectors(:,2),scaling,'r');
            qh3 = quiver(xspace(:),yspace(:),AA_shear_vectors(:,1),AA_shear_vectors(:,2),scaling,'b');
%             qh4 = quiver(xspace(:),yspace(:),-AA_shear_vectors(:,1),-AA_shear_vectors(:,2),scaling,'b');
            set(qh1,'LineWidth',1);
%             set(qh2,'LineWidth',1);
            set(qh3,'LineWidth',1);
%             set(qh4,'LineWidth',1);
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Principal simple shear quiver 4: Superimposed shear axes');
            axis equal
            
             % (5) Both orientations, reduced
             downsample_factor = 5;
             reduced_scaling = 0.8;  % These are good settings for DS18
%             downsample_factor = 2;
%             reduced_scaling = 0.5;  % These are good settings for DS2
            figure;
            xspace_reduced = xspace(1:downsample_factor:end,1:downsample_factor:end);
            yspace_reduced = yspace(1:downsample_factor:end,1:downsample_factor:end);
            ABv1 = reshape(AB_shear_vectors(:,1),size(xspace));
            ABv2 = reshape(AB_shear_vectors(:,2),size(xspace));
            AAv1 = reshape(AA_shear_vectors(:,1),size(xspace));
            AAv2 = reshape(AA_shear_vectors(:,2),size(xspace));
            ABv1_red = ABv1(1:downsample_factor:end,1:downsample_factor:end);
            ABv2_red = ABv2(1:downsample_factor:end,1:downsample_factor:end);
            AAv1_red = AAv1(1:downsample_factor:end,1:downsample_factor:end);
            AAv2_red = AAv2(1:downsample_factor:end,1:downsample_factor:end);
            AB_shear_vectors_reduced = [ABv1_red(:),ABv2_red(:)];
            AA_shear_vectors_reduced = [AAv1_red(:),AAv2_red(:)];
            qh1 = quiver(xspace_reduced(:),yspace_reduced(:),AB_shear_vectors_reduced(:,1),AB_shear_vectors_reduced(:,2),reduced_scaling,'r');
            hold on
%             qh2 = quiver(xspace_reduced(:),yspace_reduced(:),-AB_shear_vectors_reduced(:,1),-AB_shear_vectors_reduced(:,2),reduced_scaling,'r');
            qh3 = quiver(xspace_reduced(:),yspace_reduced(:),AA_shear_vectors_reduced(:,1),AA_shear_vectors_reduced(:,2),reduced_scaling,'b');
%             qh4 = quiver(xspace_reduced(:),yspace_reduced(:),-AA_shear_vectors_reduced(:,1),-AA_shear_vectors_reduced(:,2),reduced_scaling,'b');
            lw = 1;
            set(qh1,'LineWidth',lw);
%             set(qh2,'LineWidth',lw);
            set(qh3,'LineWidth',lw);
%             set(qh4,'LineWidth',lw);
            set(qh1,'MaxHeadSize',1);
%             set(qh2,'MaxHeadSize',1);
            set(qh3,'MaxHeadSize',1);
%             set(qh4,'MaxHeadSize',1);
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Principal simple shear quiver 5: Superimposed shear axes, reduced');
            axis equal
            
        end
        
        
        
        
        
        
        function makeQuiverPlotFromBlinking(obj,saveplotflag)
            figh = figure;
            DSCx = obj.DSC_fit_storage(:,:,1);
            DSCy = obj.DSC_fit_storage(:,:,2);
            [xspace,yspace] = meshgrid(obj.xaxis,obj.yaxis);
            if numel(xspace) ~= numel(DSCx)
                [xspace,yspace] = meshgrid(1:size(DSCx,1),1:size(DSCx,2));
            end
            quiver(xspace(:),yspace(:),DSCx(:),DSCy(:),0,'k');
            xlabel('x (nm)');
            ylabel('y (nm)');
            title('Displacement field');
            %             set(gca,'ydir','reverse')
            axis equal
            
            if saveplotflag
                currentd = pwd;
                cd(obj.saved_plots_folderpath)
                if ~exist(obj.saved_plots_foldername,'dir')
                    mkdir(obj.saved_plots_foldername);
                end
                cd(obj.saved_plots_foldername);
                saveas(figh,'displacement_quiverplot.png');
                cd(currentd);
            end
        end
        
        
        
        
        function makeInteractiveDisplacementDarkFieldPlot(obj,use_extended_zone,filteramp,filterrange)
            if nargin < 3
                filteramp = [];
                filterrange = [];
            end
            if use_extended_zone
                DSC1 = obj.annealed_dfield(:,:,1);
                DSC2 = obj.annealed_dfield(:,:,2);
            else
                if isempty(filteramp)
                    DSC1 = obj.DSC_fit_storage(:,:,1);
                    DSC2 = obj.DSC_fit_storage(:,:,2);
                else
                    filtered_displacement = displacementContinuityFilter(obj,filteramp,filterrange);  % The old call.
                    DSC1 = filtered_displacement(:,:,1);
                    DSC2 = filtered_displacement(:,:,2);
                end
            end
            f = figure;
            
            set(gcf, 'Position',  [0, 0, 1200, 400])
            s1 = subplot(1,3,1);
            %             ax1 = axes('Parent',s1);
            h1 = scatter(DSC1(:),DSC2(:),3,'Filled');
            title('Interactive Displacement Plot');
            % axis tight
            ylim([-0.1,1.5]);
            xlim([-1.5,1.5]);
            xlabel('X displacement');
            ylabel('Y displacement');
            axh = gca;
            plotFullDisplacementHexagons(axh);
            
            s2 = subplot(1,3,2);
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
            %             set(s1,'position',[0.1,0.1,0.4,0.8]);
            
            s3 = subplot(1,3,3);
            
            
            e = imrect(s1);
            fcn = @(pos) obj.darkFieldCallback(pos,s2,s3);
            fcn(e.getPosition());
            id = addNewPositionCallback(e,fcn);
        end
        
        
        
        function makeInteractiveInverseDarkFieldPlot(obj,use_extended_zone,continuity_threshold,filter_range)
            if nargin < 3
                continuity_threshold = [];
                filter_range = [];
            end
            if use_extended_zone
                DSC1 = obj.annealed_dfield(:,:,1);
                DSC2 = obj.annealed_dfield(:,:,2);
            else
                if isempty(continuity_threshold)
                    DSC1 = obj.DSC_fit_storage(:,:,1);
                    DSC2 = obj.DSC_fit_storage(:,:,2);
                else
                    filtered_displacement = displacementContinuityFilter(obj,continuity_threshold,filter_range);  % The old call.
                    DSC1 = filtered_displacement(:,:,1);
                    DSC2 = filtered_displacement(:,:,2);
                end
            end
            f = figure;
            set(gcf, 'Position',  [0, 0, 1000, 400])
            s1 = subplot(1,2,1);
            %             ax1 = axes('Parent',s1);
            % These are the references so that we can see the
            h1 = scatter(DSC1(:),DSC2(:),3,'Filled','b');
            title('Interactive Displacement Plot');
            % axis tight
            if ~use_extended_zone
                ylim([-0.1,1.5]);
                xlim([-1.5,1.5]);
            else
                axis equal
            end
            xlabel('X displacement');
            ylabel('Y displacement');
            axh = gca;
            plotFullDisplacementHexagons(axh);
            if ~use_extended_zone
                [ legend_handles ] = makeConcentricRings( axh );
                lh = legend(legend_handles,'Amplitude = 0.25','Amplitude = 0.5','Amplitude = 0.75','Amplitude = 1','Amplitude = 1.25','Amplitude = 1.5');
                set(lh,'Position',[0.344829545454545 0.15 0.110000000000000 0.181250000000000]);
            end
            % NPK modification 04/08: Even if we are looking at the
            % extended zone displacement in the scatterplot, still we need
            % to visualize the reduced zone for intelligibility given our
            % current plotting functions.
            if ~isempty(continuity_threshold)
                displacement_to_plot = obj.displacementContinuityFilter(continuity_threshold,filter_range);
            else
                displacement_to_plot = obj.DSC_fit_storage;
            end
            OLD = false;
            SP_unique_flag = 1;
            [ RGB_color_stack ] = getCustomDisplacementColor( displacement_to_plot, [], [], SP_unique_flag );
            s2 = subplot(1,2,2);
            %             ax2 = axes('Parent',s2);
            if OLD
                ampfield = sqrt(obj.DSC_fit_storage(:,:,1).^2 + obj.DSC_fit_storage(:,:,2).^2);
                colormap(fire);
                % Temporary fix for compatibility with partial fit
                % % %             h2 = imagesc(obj.xaxis,obj.yaxis,ampfield);
                h2 = imagesc(ampfield);
                colorbar;
            else
                h2 = imagesc(RGB_color_stack);
            end
            
            set(s2,'ydir','normal');
            axis equal
            % % %             xlim([obj.xaxis(1)-0.5*obj.scan_stepsize, obj.xaxis(end)+0.5*obj.scan_stepsize]);
            % % %             ylim([obj.yaxis(1)-0.5*obj.scan_stepsize, obj.yaxis(end)+0.5*obj.scan_stepsize]);
            xlabel('x (nm)');
            ylabel('y (nm)');
            %             set(gca,'position',[0.02 0.04 0.96 0.92]);
            set(s1,'position',[0.1,0.1,0.4,0.8]);
            
            % In the inverse plot, the roi rectangle goes on the
            % amplitude plot.
            e = imrect(s2);
            fcn = @(pos) obj.inverseDarkFieldCallback(pos,s1,use_extended_zone);
            fcn(e.getPosition());
            id = addNewPositionCallback(e,fcn);
            
            %             uib = uibutton('Text','Add a new rectangle.');
            %             fcn2 = @(pos) obj.addRect();
            %             pb = uicontrol('Style','pushbutton','Callback',@obj.pushbutton1_Callback);
            
        end
        
        
        
        
        % Do not call directly with roi_data. This is for the
        % visualizeROIIntegration function only.
        function visualizeIntegration(obj,roi_data)
            if nargin < 2
                [rs,cs,hs] = size(obj.disk_averages);
                data_to_use = reshape(obj.disk_averages,[rs*cs,1,hs]);
                data_to_use = squeeze(data_to_use);
            else
                data_to_use = double(roi_data);
            end
            f = figure;
            r = 2;
            c = 6;
            maxval = max(max(data_to_use));
            minval = min(min(data_to_use));
            edges = linspace(1.05*minval,1.05*maxval,51);
            for i = 1:12
                subplot(r,c,i);
                h = histogram(data_to_use(:,i),edges);
                h.EdgeColor = h.FaceColor;
                title(sprintf('Disk %d',i));
                xlabel('Average Pixel Counts in Disk');
                ylabel('# of disks');
                hold on
            end
            set(f,'Position',[0,600,1300,300]);
            for i = 1:12
                subplot(r,c,i);
                yvals = get(gca,'YLim');
                this_prefactor = obj.trig_prefactors(i);
                line([this_prefactor,this_prefactor],[yvals(1),yvals(2)],'Color','r');
                ylim(yvals);
            end
            
            f2 = figure;
            r = 2;
            c = 6;
            edges2 = linspace(-1,1,51);
            for i = 1:12
                subplot(r,c,i);
                h = histogram(data_to_use(:,i)/maxval,edges2);
                h.EdgeColor = h.FaceColor;
                title(sprintf('Disk %d globally normalized',i));
                xlabel('Normalized Pixel Counts');
                ylabel('# of disks');
            end
            for i = 1:12
                subplot(r,c,i);
                yvals = get(gca,'YLim');
                this_prefactor = obj.trig_prefactors(i);
                line([1,1],[yvals(1),yvals(2)],'Color','r');
                ylim(yvals);
                xlim([0,1.05]);
            end
            set(f2,'Position',[0,300,1300,300]);
            
            f3 = figure;
            r = 2;
            c = 6;
            for i = 1:12
                subplot(r,c,i);
                thismax = max(data_to_use(:,i));
                h = histogram(data_to_use(:,i)/thismax,edges2);
                h.EdgeColor = h.FaceColor;
                title(sprintf('Disk %d diskwise normalized',i));
                xlabel('Normalized Pixel Counts');
                ylabel('# of disks');
            end
            for i = 1:12
                subplot(r,c,i);
                yvals = get(gca,'YLim');
                this_prefactor = obj.trig_prefactors(i);
                line([1,1],[yvals(1),yvals(2)],'Color','r');
                ylim(yvals);
                xlim([0,1.05]);
            end
            set(f3,'Position',[0,0,1300,300]);
            
            % NPK 03/25/20: normalize based on trig height
            f4 = figure;
            r = 2;
            c = 6;
            edges3 = linspace(-1,2,51);
            for i = 1:12
                subplot(r,c,i);
                h = histogram(data_to_use(:,i)/obj.trig_prefactors(i),edges3);
                h.EdgeColor = h.FaceColor;
                title(sprintf('Disk %d',i));
                xlabel('Normalized to prefactor');
                ylabel('# of disks');
            end
            for i = 1:12
                subplot(r,c,i);
                yvals = get(gca,'YLim');
                this_prefactor = obj.trig_prefactors(i);
                line([1,1],[yvals(1),yvals(2)],'Color','r');
                ylim(yvals);
            end
            set(f3,'Position',[0,0,1300,300]);
        end
        
        
        
        
        function visualizeROIIntegration(obj,text)
            %             figh = obj.makeBlinkingPlots([],0,'flat');
            %             numchoice = input('Please enter the number of the disk you would like to use for defining AB stacking.');
            %             this_disk = obj.disk_averages(:,:,numchoice);
            if strcmp(text,'all')
                [r,c,h] = size(obj.disk_averages);
                roi_data = zeros(r*c,h);
                for i = 1:12
                    disk_data = obj.disk_averages(:,:,i);
                    roi_data(:,i) = disk_data(:);
                end
            else
                [ roi_data ] = obj.roiBoilerplate();
            end
            obj.visualizeIntegration(roi_data);
        end
        
        
        
        
        function visualizeROIDetectorValues(obj)
            % Make mask in reciprocal space, not real space.
            idx1 = 2;
            idx2 = 2;
            obj.plotSingleDP(idx1,idx2,0,1,1,1,0);
            close
            f = gcf;
            reciprocalBW = roipoly();
            
            loadnums = [1,obj.num_load_chunks-1];
            data1 = obj.partialLoad(loadnums(1));
            data2 = obj.partialLoad(loadnums(2));
            data = cat(3,data1,data2);
            [r,c,h1,h2] = size(data);
            
            histogram_values = zeros(1,0);
            for i = 1:h1
                for j = 1:h2
                    thisDP = data(:,:,i,j);
                    %                     this_histvals = thisDP(reciprocalBW & ~(obj.hBN_mask | obj.beamstop_mask));%obj.graphene_mask | obj.beamstop_mask));
                    this_histvals = thisDP(reciprocalBW & ~(obj.hBN_mask | obj.graphene_mask | obj.beamstop_mask));
                    histogram_values = vertcat(histogram_values,this_histvals);
                end
            end
            
            figure;
            histogram(histogram_values)
            h = histogram(histogram_values);
            meanval = mean(double(histogram_values))
            
            %             [ log_lik ] = PoissonFittingFunction( lambda,values )
            options = optimset('fminsearch');
            options.Display = 'iter';
            poisfun = @(lambda) PoissonFittingFunction( lambda,double(histogram_values) );
            pois_best_lambda = fminsearch(poisfun,meanval,options)
            xmax = get(gca,'XLim');
            xmax = xmax(2);
            support = 0:xmax;
            Pois_vals = poisspdf(support,pois_best_lambda);
            maxpois = max(Pois_vals);
            maxdata = max(h.Values);
            new_Pois_vals = Pois_vals.*double(maxdata)./maxpois;
            hold on;
            plot(support,new_Pois_vals,'ro-');
            
        end
        
        
        
        
        
        function [amplitude_data,binned_amplitude_counts,bin_centers] = makeDisplacementHistogram(obj,filteramp,filterrange)
            if isempty(filteramp)
                filtered_displacement = obj.DSC_fit_storage;
            else
                filtered_displacement = obj.displacementContinuityFilter(filteramp,filterrange);
            end
            amplitude = sum(filtered_displacement.^2,3).^0.5;
            figure
            histogram(amplitude(:));
            xlabel('Displacement amplitude (Angstrom)');
            ylabel('Counts');
            set(gca,'FontSize',14);
            amplitude_data = amplitude(:);
            h = histogram(amplitude(:));
            bin_centers = movmean(h.BinEdges,2);
            bin_centers(1) = [];  % this is just the first element.
            binned_amplitude_counts = h.Values;
        end
        
        
        
        
        
        function trigheights = setTrigHeight(obj)
            [ roi_data ] = obj.roiBoilerplate();
            meanvals = mean(roi_data);
            meanvals(isnan(meanvals)) = 0;
            % Correct height for the trig functions under the assumption
            % of a truly
            trigheights = meanvals .* [4 4 4 4 4 4 1 1 1 1 1 1];
            obj.trig_prefactors = trigheights;
        end
        
        
        
        
        % Saves the current instance of the FourDSTEM_Analysis_Engine to
        % the results folder so that you do not have to re-run the script
        % to manipulate the fitting results.
        function saveObject(obj)
            name = inputname(1);
            command = sprintf('save(''objectdata.mat'',''%s'')',name);
            
            currentd = pwd;
            cd(obj.saved_plots_folderpath)
            if ~exist(obj.saved_plots_foldername,'dir')
                mkdir(obj.saved_plots_foldername);
            end
            cd(obj.saved_plots_foldername);
            evalin('base',command);
            cd(currentd);
        end
        
        
        
        
        function [all_disk_mask,inner_disk_mask,outer_disk_mask,mask_stack] = getMasksForDiskIntegration(obj)
            inner_disk_mask = false(obj.datacube_size(1),obj.datacube_size(2));
            outer_disk_mask = false(obj.datacube_size(1),obj.datacube_size(2));
            mask_stack = false(obj.datacube_size(1),obj.datacube_size(2),12);
            xbase = 1:obj.datacube_size(1);
            ybase = 1:obj.datacube_size(2);
            [xspace,yspace] = meshgrid(xbase,ybase);
            for q = 1:6
                x0 = obj.disk_centers(q,1);
                y0 = obj.disk_centers(q,2);
                this_imask = isInCircle(xspace,yspace,x0,y0,obj.disk_radius);
                inner_disk_mask = inner_disk_mask | this_imask;
                mask_stack(:,:,q) = this_imask;
            end
            for q = 7:12
                x0 = obj.disk_centers(q,1);
                y0 = obj.disk_centers(q,2);
                this_omask = isInCircle(xspace,yspace,x0,y0,obj.disk_radius);
                outer_disk_mask = outer_disk_mask | this_omask;
                mask_stack(:,:,q) = this_omask;
            end
            all_disk_mask = inner_disk_mask | outer_disk_mask;
        end
        
        
        
        
        
        function makeMasksForElasticFit(obj)
            % Modified by NPK on 03/12/2020 to take advantage of the
            % smoother image obtained by averaging over the dataset.
            %             [data] = double(obj.singleLoad(2,2));  % Sometimes the first one is unusually bright for some reason.
            data = obj.averaged_DP;
            % For the hBN
            if isempty(obj.elastic_fit_mask)
                
                disp('Prepare to mask off the hBN lattice and the beamstop.');
                [ ~, hBN_mask, ~, beamstop_mask, ~ ] = makeHexagonalLatticeMask( data );
                obj.hBN_mask = hBN_mask;
                obj.beamstop_mask = beamstop_mask;
                disp('Prepare to mask off the graphene lattice.');
                [ ~, graphene_mask, ~, ~, ~ ] = makeHexagonalLatticeMask( {data} );
                obj.graphene_mask = graphene_mask;
                tf = input('Does the detector have a field of view to mask?');
                if tf
                    disp('Prepare to demarcate the detector field of view');
                    disp('Click once on the detector center and once for the detector radius.');
                    [xc,yc] = ginput(2);
                    imagesize = size(data);
                    rowbase = 1:imagesize(1);
                    colbase = 1:imagesize(2);
                    [rowspace,colspace] = meshgrid(rowbase,colbase);
                    r = sqrt((xc(1)-xc(2))^2 + (yc(1)-yc(2))^2);
                    detector_mask = ~isInCircle(colspace,rowspace,yc(1),xc(1),r);
                else
                    detector_mask = false(size(obj.graphene_mask));
                end
                obj.detector_mask = detector_mask;
                obj.elastic_fit_mask = hBN_mask | beamstop_mask | graphene_mask | detector_mask;
                show_underlying = 1;
                obj.plotElasticMask(show_underlying);
            end
        end
        
        
        
        
        function setAveragedDP(obj,loadnums)
            
            if nargin < 2
                loadnums = 1;
            end
            if strcmp(loadnums,'all')
                loadnums = 1:obj.num_load_chunks;
            elseif strcmp(loadnums,'group1')
                loadnums = 1:10:obj.num_load_chunks;
            elseif strcmp(loadnums,'group2')
                loadnums = 5:10:obj.num_load_chunks;
            end
            storedDPs = zeros(obj.datacube_size(1),obj.datacube_size(1),numel(loadnums));
            for i = 1:numel(loadnums)
                alldata = double(obj.partialLoad(loadnums(i)));
                storedDPs(:,:,i) = mean(mean(alldata,4),3);
            end
            obj.averaged_DP = mean(storedDPs,3);
        end
        
        
        
        
        
        function [com_coords] = getCOMbeamCenter(obj)
            if isempty(obj.skipped_disks)
                diskc = obj.disk_centers;
            else
                %                 error('Need to program in logic for COM with skipped disks.');
                group1 = obj.skipped_disks(obj.skipped_disks <= 6);
                group2 = obj.skipped_disks(obj.skipped_disks > 6);
                group1_res = group1-3;
                group2_res = group2-3;
                group1_res(group1_res < 1) = group1_res(group1_res < 1) + 6;
                group2_res(group2_res < 7) = group1_res(group1_res < 7) + 6;
                to_ignore = horzcat(group1,group2,group1_res,group2_res);
                diskc = obj.disk_centers;
                diskc(to_ignore,:) = [];
            end
            com_coords = mean(diskc);  % gets x,y along first dimension averaging.
            obj.com_coords = com_coords;
        end
        
        
        
        
        % Allowable function_type_strings are as follows:
        %   'power law'
        %   'lorentzian'
        %   'gaussian'
        %   'pseudo-voigt'
        %   'pseudo-voigt unconstrained'
        %   'lorentz + sinc'
        %   'ringed gaussian'
        function fitElasticScattering(obj,function_type_string)
            if isempty(obj.averaged_DP)
                obj.setAveragedDP();
            end
            data = obj.averaged_DP;
            data(obj.elastic_fit_mask) = -1;
            
            
            % Build initial parameter information on the basis of the
            % function chosen.
            y0_init = round(mean([1,size(data,1)]));
            x0_init = round(mean([1,size(data,2)]));
            use_pswarm = 0;
            const_vals = [];
            switch function_type_string
                case 'power law'
                    log10A_init = 15;
                    B_init = -4;
                    c_initial_guesses = [x0_init,y0_init,log10A_init,B_init];
                case 'lorentzian'
                    A_init = 3;
                    B_init = 50;
                    c_initial_guesses = [x0_init,y0_init,A_init,B_init];
                case 'gaussian'
                    A_init = 300;
                    B_init = 20;
                    c_initial_guesses = [x0_init,y0_init,A_init,B_init];
                case 'pseudo-voigt'
                case 'pseudo-voigt unconstrained'
                    x0_init = 206;
                    y0_init = 187;
                    AG_init = 440.6286/2;
                    AL_init = 671.5774/2;
                    fG_init = 144.3997;
                    fL_init = 116.0446;
                    c_initial_guesses = [x0_init,y0_init,AG_init,AL_init,fG_init,fL_init];
                    use_pswarm = 1;
                case 'lorentz + sinc'
                case 'ringed gaussian'
                    [com_coords] = obj.getCOMbeamCenter();
                    const_vals = com_coords;
                    %                     I0_init = 0.01;
                    %                     Ipeak_init = 10;
                    %                     Iring_init = 1;
                    %                     rRing_init = 100;
                    %                     fwhm_peak_init = 100;
                    %                     fwhm_ringinner_init = 100;
                    %                     fwhm_ringouter_init = 100;
                    I0_init = 0.01;
                    Ipeak_init = 10;
                    Iring_init = 2;
                    rRing_init = 60;
                    fwhm_peak_init = 50;
                    fwhm_ringinner_init = 20;
                    fwhm_ringouter_init = 20;
                    
                    c_initial_guesses = [I0_init,Ipeak_init,Iring_init,rRing_init,...
                        fwhm_peak_init,fwhm_ringinner_init,fwhm_ringouter_init];
            end
            
            use_fminsearch = 0;
            if use_fminsearch
                elasticFit = @(c) getFullRadialElasticFit(c,data,function_type_string,const_vals);
            else
                elasticFit = @(c) getFullRadialElasticFit_lsqnonlin(c,data,function_type_string,const_vals);
            end
            options = optimset;
            options.Display = 'iter';
            options.TolFun = 1e-7;
            options.TolX = 1e-7;
            options.MaxFunEvals = 5000;
            options.MaxIter = 5000;
            
            if use_fminsearch
                c_optimal = fminsearch(elasticFit,c_initial_guesses,options);
            else
                c_optimal = lsqnonlin(elasticFit,c_initial_guesses,[-0.1,0,0,0,0,0,0],[inf,inf,inf,inf,inf,inf,inf],options);
            end
            obj.elastic_fit_parameters = c_optimal;
            optimal_fit = getFullRadialElasticPred(c_optimal,size(data),function_type_string,const_vals);
            obj.optimal_elastic_fit = optimal_fit;
            
            figh = figure; subplot(2,2,1);
            nan_masked_data = data;
            nan_masked_data(data == -1) = nan;
            surf(nan_masked_data); shading flat;
            title('Data to fit');
            xlabel('x');
            ylabel('y');
            colormap(pink);
            
            subplot(2,2,2);
            residuals = data - optimal_fit;
            
            obj.elastic_residuals = residuals;
            residuals(data == -1) = nan;
            residuals_for_RMS = residuals(:);
            residuals_for_RMS(isnan(residuals_for_RMS)) = [];
            RMSresiduals = rms(rms(residuals_for_RMS))
            surf(residuals); shading flat;
            colormap(pink);
            title('Fitting residuals');
            xlabel('x');
            ylabel('y');
            
            
            subplot(2,2,4);
            surf(log10(abs(optimal_fit))); shading flat;
            colormap(pink);
            title(sprintf('log10 of optimal %s fit',function_type_string));
            xlabel('x');
            ylabel('y');
            
            subplot(2,2,3);
            surf(optimal_fit); shading flat;
            colormap(pink);
            title(sprintf('Optimal %s fit',function_type_string));
            xlabel('x');
            ylabel('y');
            
            set(figh,'Position',[0,0,1200,800]);
        end
        
        
        
        
        
        % Acceptible values for center_type are "residuals fit" or "com"
        function setLineAverageElasticResiduals(obj,center_type)
            switch center_type
                case 'original fit'
                    center = obj.elastic_fit_parameters(1:2);
                case 'residuals fit'
                    center = obj.elastic_residuals_fit_parameters(1:2);
                case 'com'
                    center = obj.com_coords;
            end
            nlines = 1000;
            [radius_values, line_averages] = obj.getLineAverageElasticResiduals(center,nlines,1);
            obj.elastic_residuals_radial = [radius_values,line_averages];
            obj.radial_correction_center = center;
        end
        
        
        
        
        function [radius_values, line_averages] = getLineAverageElasticResiduals(obj,center,nlines,plotflag)
            %             disp('getLineAverageElasticResiduals()');
            % Prepare the residuals matrix for interpolation
            [r,c] = size(obj.averaged_DP);
            %             center = obj.elastic_fit_parameters(1:2);
            residuals = obj.averaged_DP - obj.optimal_elastic_fit;
            residuals(obj.elastic_fit_mask) = nan;
            xbase = (1:r) - center(1);
            ybase = (1:c) - center(2);
            [xspace,yspace] = meshgrid(xbase,ybase);
            
            % Evaluate interpolation for radial lines
            %             nlines = 1000;
            angles = linspace(0,2*pi,nlines);
            nradii = 2*max(r,c);
            radius_values = linspace(0,max(r,c),nradii);
            line_storage = zeros(nradii,3,nlines);
            for i = 1:nlines
                this_angle = angles(i);
                Xq = radius_values.*cos(this_angle);
                Yq = radius_values.*sin(this_angle);
                Vq = interp2(xspace,yspace,residuals,Xq,Yq,'nearest');
                line_storage(:,:,i) = [Xq',Yq',Vq'];
            end
            nonnanvals = sum(~isnan(line_storage(:,3,:)),3);
            line_storage(isnan(line_storage)) = 0;
            line_storage_to_sum = line_storage(:,3,:);
            sum_line = sum(line_storage_to_sum,3);
            sum_line(nonnanvals == 0) = 0;  % no information there with which to correct the diffraction pattern.
            av_line = zeros(size(sum_line));
            av_line(nonnanvals > 0) = sum_line(nonnanvals > 0)./nonnanvals(nonnanvals > 0);
            radius_values = radius_values';  % for consistent return dimensions
            line_averages = av_line;
            %             obj.elastic_residuals_radial = [radius_values',av_line(:,3)];
            
            if plotflag
                figure
                imagesc(residuals);
                hold on;
                colormap('jet');
                v = viscircles(repmat(center,3,1),[50;75;100],'Color','k');
                legend([v],'Radii = 50, 75, 100');
                colorbar
                axis equal
                title('Residuals before radial subtraction');
                
                figure
                plot(radius_values,av_line);
                title('Radial fit residuals average.');
                xlabel('Radius from fitted beam center');
                ylabel('Mean magnitude of residuals');
            end
        end
        
        
        
        
        
        % Goal of this function is to find the center of the residuals ring
        % that is characteristically left over after the beamstop
        % scattering. Can conceivably do this either by maximizing the
        % absolute value of the averages (indicating that systematic
        % variations are concentrated onto particular locations in the
        % radii) or by minimizing the norm of the residuals matrix after
        % this is subtracted back off. At first, we will take the former
        % approach.
        function fitElasticResidualsCenter(obj,obj_fcn_choice)
            if obj_fcn_choice == 1
                obj_fcn = @(center_guess) obj.fitElasticResidualsCenter_AbsHelperFcn(center_guess);
            end
            options = optimset('fminsearch');
            options.MaxFunEvals = 100;
            options.Display = 'iter';
            init_center_guess = obj.elastic_fit_parameters(1:2)+[30,5];  % from the peak fit, which is usually close but not great.
            residuals_center = fminsearch(obj_fcn,init_center_guess,options);
            obj.elastic_residuals_fit_parameters = residuals_center;
        end
        
        function merit_val = fitElasticResidualsCenter_AbsHelperFcn(obj,center_guess)
            nlines = 100;
            plotflag = 0;
            [~, line_averages] = obj.getLineAverageElasticResiduals(center_guess,nlines,plotflag);
            merit_val = -norm(abs(line_averages));  % Because we want to maximize this quantity, take the negative of it.
        end
        
        
        
        % This is an interpolation through the elastic fit residuals. It
        % interpolates on both cartesian and polar coordinates
        function interpolateThroughPeakResiduals(obj,use_post_radial_fit)
            % Cartesian interpolation (suboptimal now that we have polar)
            if use_post_radial_fit
                elastic_to_subtract = obj.correctDP(obj.elastic_residuals,0,...
                    1,0,0);
                %              correctDP(obj,DP,subtract_elastic_flag,...
                %                 radial_elastic_correction_flag,cartesian_interpolation_flag,polar_interpolation_flag)
            else
                elastic_to_subtract = obj.elastic_residuals;
            end
            
            xbase = 1:obj.datacube_size(1);
            ybase = 1:obj.datacube_size(2);
            [xspace_org,yspace_org] = meshgrid(xbase,ybase);
            xspace = xspace_org;
            yspace = yspace_org;
            elastic_resid = elastic_to_subtract;
            elastic_resid(obj.elastic_fit_mask) = [];
            xspace(obj.elastic_fit_mask) = [];
            yspace(obj.elastic_fit_mask) = [];
            sI = scatteredInterpolant(xspace',yspace',elastic_resid','natural');
            interpVals = sI(xspace_org,yspace_org);
            figure;
            imagesc(interpVals);
            colormap(jet);
            colorbar;
            shading flat;
            axis square
            obj.interpolated_residuals = interpVals;
            
            % Polar interpolation
            xspace_centered = xspace_org - obj.radial_correction_center(1);
            yspace_centered = yspace_org - obj.radial_correction_center(2);
            [theta,rho] = cart2pol(xspace_centered,yspace_centered);
            elastic_resid = elastic_to_subtract;
            tref = vertcat(theta(:)+2*pi,theta(:),theta(:)-2*pi);
            rref = repmat(rho(:),3,1);
            vref = repmat(elastic_resid(:),3,1);
            reduced_mask = obj.hBN_mask(:) | obj.graphene_mask(:) | obj.beamstop_mask(:);
            mask_remove = repmat(reduced_mask,3,1);
            tquery = theta(:);
            tquery(~reduced_mask) = [];
            rquery = rho(:);
            rquery(~reduced_mask) = [];
            tref(mask_remove) = [];
            rref(mask_remove) = [];
            vref(mask_remove) = [];
            F = scatteredInterpolant(tref,rref,vref,'linear');
            interpolated_query = F(tquery,rquery);
            % get coords in the original image based on linear indexing
            % because the array sizes have remained constant.
            obj.interpolated_polar_residuals = elastic_resid;
            obj.interpolated_polar_residuals(reduced_mask) = interpolated_query;
            % optional? Doing this will remove the interpolation under the
            % beamstop, but that tends to go wonky anyways on the polar
            % interpolation
            obj.interpolated_polar_residuals(obj.beamstop_mask) = 0;
            figure;
            imagesc(obj.interpolated_polar_residuals);
            colormap(jet);
            colorbar;
            shading flat;
            axis square
        end
        
        
        
        
        function plotElasticResiduals(obj)
            figure;
            to_plot = obj.elastic_residuals;
            to_plot(obj.elastic_fit_mask) = 0;
            pcolor(to_plot);
            shading flat;
            axis square;
            colormap jet
            colorbar
            hold on;
            if numel(obj.elastic_fit_parameters) == 7  % then we used the fixed ringed Gaussian fit
                s1 = scatter(mean(obj.disk_centers(:,1)),mean(obj.disk_centers(:,2)),'k','filled');
            else
                s1 = scatter(obj.elastic_fit_parameters(1),obj.elastic_fit_parameters(2),'k','filled');
            end
            if ~isempty(obj.elastic_residuals_fit_parameters)
                s2 = scatter(obj.elastic_residuals_fit_parameters(1),obj.elastic_residuals_fit_parameters(2),'r','filled');
                legend([s1,s2],'Peak fit center','Residuals fit center');
            else
                legend([s1],'Peak fit center');
            end
            
        end
        %         thisDP = obj.correctDP(thisDP,subtract_elastic_flag,radial_elastic_correction_flag);
        % if radial_elastic_correction_flag
        %                                 % compute radius based on indices
        %                                 radius = sqrt((x0 - center(1)).^2 + (y0 - center(2)).^2);
        %                                 [~,idx] = min(abs(obj.elastic_residuals_radial(:,1)-radius));
        %                                 correction_val = obj.elastic_residuals_radial(idx,2);
        %                             end
        function correctedDP = correctDP(obj,DP,subtract_elastic_flag,...
                radial_elastic_correction_flag,cartesian_interpolation_flag,polar_interpolation_flag)
            thisDP = double(DP);
            if subtract_elastic_flag
                correctedDP = thisDP - obj.optimal_elastic_fit;
            else
                correctedDP = thisDP;
            end
            if radial_elastic_correction_flag
                METHOD = 2;
                center = obj.radial_correction_center;
                xbase = (1:obj.datacube_size(2)) - center(1);
                ybase = (1:obj.datacube_size(1)) - center(2);
                [xspace,yspace] = meshgrid(xbase,ybase);
                radii = sqrt(xspace.^2 + yspace.^2);
                if METHOD == 1
                    for i = 1:obj.datacube_size(1)
                        for j = 1:1:obj.datacube_size(2)
                            radius = radii(i,j);
                            [~,idx] = min(abs(obj.elastic_residuals_radial(:,1)-radius));
                            correction_val = obj.elastic_residuals_radial(idx,2);
                            correctedDP(i,j) = correctedDP(i,j) - correction_val;
                        end
                    end
                elseif METHOD == 2
                    % Equivalent results, but this is waaaay faster
                    vq = interp1(obj.elastic_residuals_radial(:,1),obj.elastic_residuals_radial(:,2),radii);
                    correctedDP = correctedDP - vq;
                end
            end
            if cartesian_interpolation_flag
                correctedDP = correctedDP - obj.interpolated_residuals;
            end
            % put the logic for this here, because we only want to correct
            % after the peak subtraction.
            if polar_interpolation_flag
                correctedDP = correctedDP - obj.interpolated_polar_residuals;
            end
            assert(~(cartesian_interpolation_flag && polar_interpolation_flag),...
                'Do not use both cartesian and polar residuals interpolation, as this double-counts the residual subtraction.');
        end
        
        
        
        
        
        % Visualize the correction function to ensure it is removing the
        % elastic scattering as predicted.
        function plotBaselinedData(obj,exponent,removeHBN,removeGraphene,removeBeamstop,...
                subtractScattering,subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation)
            if nargin < 2
                exponent = 1;
            end
            %             subtract_elastic_flag = 1;
            %             radial_elastic_correction_flag = 0;
            %             interpolation_correction_flag = 1;
            correctedDP = obj.correctDP(obj.averaged_DP,subtractScattering,...
                subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation);
            if removeHBN
                correctedDP(obj.hBN_mask) = 0;
            end
            if removeGraphene
                correctedDP(obj.graphene_mask) = 0;
            end
            if removeBeamstop
                correctedDP(obj.beamstop_mask) = 0;
            end
            figure
            surf(correctedDP); shading flat;
            colormap(gray)
            
            plotDP(correctedDP,exponent);
            colormap(gray);
            hold on;
%             scatter(obj.elastic_fit_parameters(1),obj.elastic_fit_parameters(2),'b','filled');
            
        end
        
        
        
        function plotRadialResiduals(obj)
            figure;
            plot(obj.elastic_residuals_radial(:,1),obj.elastic_residuals_radial(:,2));
            xlabel('Radius from beam center');
            ylabel('Magnitude of averaged residuals');
            title(sprintf('Beam center is %d, %d',obj.elastic_fit_parameters(1),obj.elastic_fit_parameters(2)));
        end
        
        
        
        
        function plotElasticMask(obj,show_underlying)
            % Plot the mask as a whole
            if isempty(obj.averaged_DP)
                [data] = obj.singleLoad(2,2);
            else
                data = obj.averaged_DP;
            end
            
            if show_underlying
                exponent = 0.3;
                factor = 10;
                filtered_image = double(data).^exponent/factor;
                filtered_image(obj.elastic_fit_mask) = filtered_image(obj.elastic_fit_mask)*factor;
                % f2 = figure;
                % imagesc(filtered_image);
            else
                filtered_image = double(data);
                filtered_image(obj.elastic_fit_mask) = 0;
                exponent = 1;
            end
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
        
        
        
        
        
        function [f,objfcn] = computeRMSRSurface(obj,indices,extended_zone_flag,individual_residuals_flag)
            displacement = permute(obj.DSC_fit_storage(indices(1),indices(2),:),[1,3,2]);
            integration_values = permute(obj.disk_averages(indices(1),indices(2),:),[3,1,2]);
            weight_vector_local = obj.weight_vector;
            objfcn = @(DSCvector) (integration_values' - trigFittingFunctions(DSCvector,obj.trig_prefactors,extended_zone_flag)).*weight_vector_local;
            if extended_zone_flag
                xbase = -4:0.05:4;
                ybase = -4:0.05:4;
            else
                xbase = -1.5:0.01:1.5;
                ybase = 0:0.01:1.5;
            end
            [xspace,yspace] = meshgrid(xbase,ybase);
            RMSR = zeros(size(xspace));
            residuals = zeros([size(xspace),12]);
            for i = 1:numel(ybase)
                if ~mod(i,10)
                    fprintf('Computing RMSR surface slice %d out of %d.\n',i,numel(ybase));
                end
                for j = 1:numel(xbase)
                    this_DSC = [xspace(i,j),yspace(i,j)];
                    result = objfcn(this_DSC);  % spits out the vector of residuals
                    computeresult = result;
                    computeresult(obj.skipped_disks) = [];
                    RMSR(i,j) = rms(computeresult);
                    residuals(i,j,:) = permute(result,[1,3,2]);
                end
            end
            
            f = figure;
            if extended_zone_flag
                contourf(xbase,ybase,RMSR,15);
            else
                contourf(xbase,ybase,RMSR,30);
            end
            title(sprintf('RMSR surface, indices [%d,%d]',indices(1),indices(2)),'FontSize',12);
            axis equal
            hold on
            if extended_zone_flag
                plotFullDisplacementHexagons(gca,3);
            else
                plotFullDisplacementHexagons(gca,1);
                xlim([-1.5, 1.5]);
                ylim([-0.1, 1.7]);
            end
            colormap jet
            h = colorbar;
            ylabel(h, 'Root-mean-square interferometry residuals');
            xlabel('x displacement (angstroms)');
            ylabel('y displacement (angstroms)');
            set(gca,'FontSize',16);
            set(h,'FontSize',16);
            set(gcf,'Position',[0,0,800,800]);
        end
        
        
        
        
        function makeMultistartMovie(obj,indices)
            [f,objfcn] = obj.computeRMSRSurface(indices,0,[]);
            options.Display = 'off';
            lb = [-1.24,0];
            ub = [1.24,1.43];
            [ best_DSC, convergence_storage, video_results ] = extendedMultistartDiskOptimization( objfcn, lb, ub, options, 1, 1 );
            % filter to make a better movie
            for i = 1:numel(video_results)
                thisres = video_results{i};
                lastres = [inf,inf];
                delmarks = [];
                for j = 1:size(thisres,1)
                    if all(abs(thisres(j,:)-lastres)<1e-10)
                        delmarks = [delmarks,j];
                    end
                    lastres = thisres(j,:);
                end
                video_results{i}(delmarks,:) = [];
            end
            
            % plot/movie
            figure(f);
            colormap(pink);
            hold on;
            for i = 1:numel(video_results)
                thisres = video_results{i};
                plot(thisres(:,1),thisres(:,2),'ro-','LineWidth',2);
                scatter(thisres(1,1),thisres(1,2),150,[0.7,0,0],'filled');
                plot(thisres(end,1),thisres(end,2),'wp','MarkerSize',12);
            end
        end
        
        
        
        
        
        % Because it looks like we'll be using this algorithm rather
        % extensively, need to have a way of storing the results for a good
        % run (maybe also a struct specifying the parameter values chosen)
        function setMultistartFilteredData(obj)
            
        end
        
        
        
        
        
        % Note that this method is always going to use the annealed
        % displacement field. The real question is whether you want to do
        % total variation smoothing first.
        function [nondimensionalized_displacement_line,TwoD_displacements] = compute1DSolitonShearStrainFromROI(obj,saveplot_flag,filterstruct,fix_zero_one_points,cut_number)
            if nargin < 5
                cut_number = [];
            end
            
            if isempty(obj.strain_1D_roi_storage_cell)
                obj.strain_1D_roi_storage_cell = cell(1,0);
                obj.strain_1D_percent_values = zeros(1,0);
                obj.soliton_width_1D_values = zeros(1,0);
            else
                disp('It appears some soliton boundary 1D strain values have already been computed.');
                yn = input('Do you want to (1) overwrite these values or (0) append to them?');
                if yn
                    obj.strain_1D_roi_storage_cell = cell(1,0);
                    obj.strain_1D_percent_values = zeros(1,0);
                    obj.soliton_width_1D_values = zeros(1,0);
                end
            end
            
            if (nargin < 2), saveplot_flag = false; end
            if (nargin < 3), TV_filter_epsilons = []; end
            if (nargin < 4), fix_zero_one_points = false; end
            
            warning('Please ensure that the scan_stepsize variable is set correctly in this class! Otherwise the strain calculation will be inaccurate.');
            obj.reformAxes();
            
            % Changed by NPK on 04/05/2020 to be the annealed
            % displacement field because that's what we will always use
            % for strain calculations.
            dfieldx = obj.annealed_dfield(:,:,1);  % because we are referencing this against the v1
            dfieldy = obj.annealed_dfield(:,:,2); % change this momentarily
            [ dfieldx,dfieldy ] = filterDisplacement( dfieldx,dfieldy,filterstruct,0,obj );
%             if ~isempty(TV_filter_epsilons)
%                 plotflag = 1;
%                 dfieldx = TVM_image_denoise( dfieldx, TV_filter_epsilons(1), plotflag, 'X displacement', parula );
%                 dfieldy = TVM_image_denoise( dfieldy, TV_filter_epsilons(2), plotflag, 'Y displacement', parula );
%                 disp('The total variation denoised displacement field will be used when computing the strain.');
%             end
            
            % Use this anonymous handle to set the displacement plotting
            % options.
            %             if isempty(TV_filter_epsilons)
%             genfig = @() obj.makeCustomDisplacementColorPlot([],[],[],[],1);
            [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( [dfieldx(:),dfieldy(:)] );
            rdx = reshape(reduced_zone_disps(:,1),size(dfieldx));
            rdy = reshape(reduced_zone_disps(:,2),size(dfieldy));
            genfig = @() obj.makeOutsideCustomDisplacementColorPlot(rdx,rdy);
            %             end
            linedrawn = false;
            terminated = false;
            count = 0;
            while ~terminated
                [fc,axc] = genfig();
                num_solitons_computed = size(obj.strain_1D_roi_storage_cell,1);
                fprintf('Preparing to calculate the %dth 1D soliton strain.\n',num_solitons_computed);
                
                disp('Please select the left and right endpoints of the SP region to analyze, parallel along the direction of the soliton wall.');
                coords = ginput(2);
                
                while true
                    r = input('Please enter the distance in each direction the fitting region may go: ');
                    if linedrawn
                        close(fc);
                        [fc,axc] = genfig();
                    end
                    figure(fc);
                    axes(axc);
                    Linep = coords(2,:) - coords(1,:);
                    parallel_angle = atan(Linep(2)/Linep(1));
                    perp_angle = pi/2 + parallel_angle;
                    hold on
                    inc = [r*cos(perp_angle),r*sin(perp_angle)];
                    bcoords1 = coords(1,:) + inc;
                    bcoords2 = coords(1,:) - inc;
                    bcoords3 = coords(2,:) + inc;
                    bcoords4 = coords(2,:) - inc;
                    
                    l1 = line([bcoords1(1),bcoords2(1)],[bcoords1(2),bcoords2(2)],'Color','g','LineWidth',2.5);
                    l2 = line([bcoords2(1),bcoords4(1)],[bcoords2(2),bcoords4(2)],'Color','g','LineWidth',2.5);
                    l3 = line([bcoords4(1),bcoords3(1)],[bcoords4(2),bcoords3(2)],'Color','g','LineWidth',2.5);
                    l4 = line([bcoords3(1),bcoords1(1)],[bcoords3(2),bcoords1(2)],'Color','g','LineWidth',2.5);
                    boundary_coords = vertcat(bcoords1,bcoords2,bcoords3,bcoords4);
                    linedrawn = true;
                    
                    yn = input('Is the displayed ROI defined satisfactorily?');
                    outercontinue = 0;
                    if ~yn
                        yn = input('Do you want to (0) modify the distance of the rectangular ROI, or (1) redefine the points of the ROI?');
                        if yn
                            outercontinue = 1;
                            break
                        else
                            outercontinue = 0;
                        end
                        
                    else
                        % save the roi figure
                        if saveplot_flag
                            fh = gcf;
                            titlename = sprintf('ROI for %dth 1D strain calculation',count);
                            title(titlename);
                            savename1 = sprintf('1DstrainROI_%d.fig',count);
                            savename2 = sprintf('1DstrainROI_%d.png',count);
                            currentd = pwd;
                            cd(obj.saved_plots_folderpath);
                            if ~exist(obj.saved_plots_foldername,'dir')
                                mkdir(obj.saved_plots_foldername);
                            end
                            cd(obj.saved_plots_foldername);
                            savefig(fh,savename1);
                            saveas(fh,savename2);
                            cd(currentd);
                        end
                        break
                    end
                end
                if outercontinue % meaning we want to go back to the definition of the ROI region
                    continue
                end
                
                %% A good roi has been computed by this point. Compute strain via numerics
                dfieldx_relative = dfieldx;
                dfieldy_relative = dfieldy;
                
                % Get a linecut average over the soliton
                linecut = @(base_t,N) [linspace(0,r*2,N)'.*cos(perp_angle),linspace(0,r*2,N)'.*sin(perp_angle)] + repmat((1-base_t)*bcoords2 + base_t*bcoords4,N,1);
                if isempty(cut_number)
                    Ncuts = 50;
                else
                    Ncuts = cut_number;
                end
                NinCut = 100;
                linecutdisps = zeros(NinCut,2,Ncuts);
                basecoords = linspace(0,1,Ncuts);
                obj.reformAxes();
                xbase = obj.xaxis;
                ybase = obj.yaxis;
                [xspace,yspace] = meshgrid(xbase,ybase);
                for i = 1:Ncuts
                    this_linecut = linecut(basecoords(i),NinCut);
                    %     dispxq = interp2(dfieldx,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
                    %     dispyq = interp2(dfieldy,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
%                     dispxq = interp2(dfieldx_relative,this_linecut(:,1),this_linecut(:,2));
%                     dispyq = interp2(dfieldy_relative,this_linecut(:,1),this_linecut(:,2));
                    dispxq = interp2(xspace,yspace,dfieldx_relative,this_linecut(:,1),this_linecut(:,2));
                    dispyq = interp2(xspace,yspace,dfieldy_relative,this_linecut(:,1),this_linecut(:,2));
                    this_line_disps = [dispxq,dispyq];
                    linecutdisps(:,:,i) = this_line_disps;
                    % save the starting and ending coordinates for later
                    % use.
                    if i == floor(Ncuts/2)
                        lower_AB_RS_point = this_linecut(1,:);
                        higher_AB_RS_point = this_linecut(end,:);
                    end
                end
                mean_linecut_disps = mean(linecutdisps,3);
                fh3 = figure;
                plot(mean_linecut_disps(:,1),mean_linecut_disps(:,2),'b-o');
                hold on
                plotFullDisplacementHexagons(gca);
                % perp_unit_vec = [cos(perp_angle),sin(perp_angle)];
                
                % Set up Sine_Gordon (modified to get rid of the symbolic
                % part)
                %                 syms x w
                %                 delu = 2/pi*atan(exp(pi*x/w));
                %                 deluf = matlabFunction(delu);
                deluf = @(w,x) 2/pi*atan(exp(pi*x./w));
                % Option here to make a plot of the soliton calculation all
                % by itself.
                %                 xs = -5:0.01:5;
                %                 w = 4*ones(size(xs));
                %                 deluvals = deluf(w,xs);
                %                 figure
                %                 plot(xs,deluvals,'ro-');
                
                % NPK TODO: here we need to introduce logic to project the
                % displacement onto the correct displacement direction
                % (i.e. which saddle point do we have?). Maybe can even use
                % some of the colorization function logic. Note we can flip
                % sign here by taking negative sign displacement vectors if
                % needed.
                
                
                % NPK added 04/06/2020: figure out the orientation and
                % location of this saddle point.
                % Convention is that we go from smaller y to larger y for
                % the fitting region, since there are no horizontal AB
                % connection lines in the image.
                mean1 = mean_linecut_disps(1:5,:);
                mean2 = mean_linecut_disps(end-4:end,:);
                if mean1(2) > mean2(2)
                    %                     t = mean1;
                    %                     mean1 = mean2;
                    %                     mean2 = t;
                    mean_linecut_disps = flipud(mean_linecut_disps);
                end
                % Get closest emitter in real space
                %                 m4.emitter_pixel_pos lower_AB_RS_point
                
                % NPK 05/20/2020 modification
                emitter_disp_pos = [obj.xaxis(obj.emitter_pixel_pos(:,1))',obj.yaxis(obj.emitter_pixel_pos(:,2))'];
                
                lowerdists = sum((emitter_disp_pos - fliplr(lower_AB_RS_point)).^2,2).^0.5;
                [~,minlidx] = min(lowerdists);
                higherdists = sum((emitter_disp_pos - fliplr(higher_AB_RS_point)).^2,2).^0.5;
                [~,minhidx] = min(higherdists);
                lower_AB_region_disp = obj.emitter_displacements(minlidx,:);
                higher_AB_region_disp = obj.emitter_displacements(minhidx,:);
                % need to get projected displacement data
                [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
                scale = hexagon_lattice_constant/sqrt(3);
                zeropoint = lower_AB_region_disp;
                zeroed_disps = mean_linecut_disps - repmat(zeropoint,NinCut,1);
                % Get vector connecting the two AB regions
                v_AB = higher_AB_region_disp - lower_AB_region_disp;
                projected_disps = v_AB./(norm(v_AB))*zeroed_disps';
                scaled_disps = projected_disps/scale;
                % if this is set, then we will use a fixed 0 - 1 range for
                % the points
                if fix_zero_one_points
                    min_disp = mean(scaled_disps(end-4:end));
                    scaled_disps = scaled_disps - min_disp;
                    max_disp = mean(scaled_disps(1:5));
                    scaled_disps = scaled_disps/max_disp;
                end
                
                
                % ascertain the distance
                if scaled_disps(end) > scaled_disps(1)
                    scaled_disps = fliplr(scaled_disps);
                end
                % 05/20/2020 Fixed this for the new case where we have
                % reformed the axes.
%                 warning('Please ensure that the scan_stepsize variable is set correctly in this class! Otherwise the strain calculation will be inaccurate.');
%                 xvals = linspace(0,2*r*obj.scan_stepsize,NinCut);
                xvals = linspace(0,2*r,NinCut);  % Now r already takes into account the stepsize.
                fdata = figure;
                plot(fliplr(xvals),scaled_disps,'bo');
                xlabel('Distance perpendicular to soliton boundary (nm)');
                ylabel('Nondimensionalized stacking change');
                
                % Conduct fit
                fitpred = @(w,translation,xs) deluf(w*ones(1,numel(xs)),xs-translation);
                fitresid = @(params) fliplr(scaled_disps) - fitpred(params(1),params(2),xvals);
                fitfun = @(params) rms(fitresid(params));
                beta0 = [r/2,r];
                options.Display = 'iter';
                options.MaxFunEvals = 1e4;
                beta = fminsearch(fitfun,beta0,options);
                % Examine results
                figure(fdata)
                predvals = fitpred(beta(1),beta(2),xvals);
                hold on
                plot(xvals,predvals,'r-');
                beta
                
                % Use the strain derivative function to get max strain across soliton boundary
                units_strain_percent = @(w_soliton) 0.07104/w_soliton * 100;
                strain_percent = units_strain_percent(beta(1));
                fprintf('The computed strain is %.3f%%.\n',strain_percent);
                soliton_width = beta(1);
                fprintf('The computed soliton width is %.3f nm.\n',soliton_width);
                
                % Store the strain calculation results and roi information
                % rows of the roi cell and the array correspond.
                obj.strain_1D_roi_storage_cell{end+1,1} = boundary_coords;  % The roi can be reconstructed. TODO refactor plotting code.
                obj.strain_1D_percent_values(end+1,1) = strain_percent;
                obj.soliton_width_1D_values(end+1,1) = beta(1);
                                
                % Make plots if desired
                if saveplot_flag
                    fh2 = gcf;
                    titlename = sprintf('%dth soliton strain fit: %.3f%%',count,strain_percent);
                    title(titlename);
                    figure(fh3);
                    titlename2 = sprintf('%dth soliton strain data to fit: %.3f%%',count,strain_percent);
                    title(titlename2);
                    savename1 = sprintf('1Dstrainfit_%d.fig',count);
                    savename2 = sprintf('1Dstrainfit_%d.png',count);
                    savename1b = sprintf('1Dstrain_avdata_%d.fig',count);
                    savename2b = sprintf('1Dstrain_avdata_%d.png',count);
                    currentd = pwd;
                    cd(obj.saved_plots_folderpath);
                    if ~exist(obj.saved_plots_foldername,'dir')
                        mkdir(obj.saved_plots_foldername);
                    end
                    cd(obj.saved_plots_foldername);
                    savefig(fh2,savename1);
                    saveas(fh2,savename2);
                    savefig(fh3,savename1b);
                    saveas(fh3,savename2b);
                    cd(currentd);
                end
                
                % Query for next steps
                count = count + 1;
                tf = input('Do you wish to compute another soliton strain value? 1/0 ');
                if ~tf
                    nondimensionalized_displacement_line = [fliplr(xvals)',scaled_disps'];
                    TwoD_displacements = mean_linecut_disps;
                    return
                end
            end
        end
        
        
        
        
        
        
        % Helper function for circlefit and ellipsefit methods for AA
        % domains 
        function [AA_centroids,amask,f] = findAAcenters(obj,filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma)
            MIN_AA_SPACING = 5;  % nm
            useNew = true;
            
            if isempty(filterrange)
                displacement_field = obj.DSC_fit_storage;
            else
                displacement_field = obj.displacementContinuityFilter(filteramp,filterrange);
            end
            if nargin < 8
                AAblur_Gaussian_sigma = 5; 
            end
            if hotspot_suppression
                % The old behavior
                if isempty(obj.hotspot_suppression_mask)
                    displacement_field = obj.suppressHotspots(displacement_field);
                else  % 05/28/2020 the new behavior targeted for dataset 22
                    displacement_field = obj.displacementContinuityFilter(filteramp,filterrange,1,[],[],hotspot_suppression);
                end
            end
            amplitude = sqrt(displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2);
            [ RGB_color_stack, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, [], [], 1, 1 );
            S = HSV_color_stack(:,:,2);
            V = HSV_color_stack(:,:,3);
            f = figure;
            imagesc(obj.xaxis,obj.yaxis,amplitude);
            title('Filtered amplitude');
            colormap(gray);
            colorbar;
            set(gca,'yDir','normal');
            [figh] = obj.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1);
            %             m4.makeCustomDisplacementColorPlot([],[],0.3,9,1);
            
            thresh = AA_amp_threshold;
            amask = amplitude < thresh;
            figure(figh);
            hold on
            h = imagesc(obj.xaxis,obj.yaxis,cat(3,~amask,~amask,zeros(size(amask))));
            h.AlphaData = 0.5;
            
            %             gauss_filter_thresh = 0.1;
            %             edge_threshold = 10;
            if useNew
                amask = bwmorph(amask,'clean');
                amask = bwareaopen(amask,5);  % clear out light spots
            end
            
            fAA = imgaussfilt(double(amask),AAblur_Gaussian_sigma);
            if useNew
                fAA(obj.tear_mask) = 0;
            end
            fAA(fAA < gauss_filter_threshold) = 0;
            
            
            
            maxvals = imregionalmax(fAA);
            lininds = find(maxvals);
            [yv,xv] = ind2sub(size(fAA),lininds);
            % Added 05/29/2020: ascertain and remove AA centers that are
            % closer together than MIN_AA_SPACING
            ydists = yv - yv';
            xdists = xv - xv';
            dists = sqrt(xdists.^2 + ydists.^2);
            flagmat = zeros(size(dists));
            flagmat(dists < MIN_AA_SPACING) = 1;
            flagmat = flagmat - eye(size(dists,1));
            if nnz(flagmat > 0)
                [row,col] = find(triu(flagmat));
                npairs = size(row,1);  % The two AA centers involved.
                to_delete = zeros(npairs,1);
                for q = 1:npairs
                    % The problem we are having is when these are exactly
                    % equal. So comparing values to decide between them
                    % isn't going to help. Just discard the second.
                    these_indices = [row(q),col(q)];
                    to_delete(q) = these_indices(2);
%                     coords1 = [yv(these_indices(1)),xv(these_indices(1))];
%                     coords2 = [yv(these_indices(2)),xv(these_indices(2))];
                end
                xv(to_delete) = [];
                yv(to_delete) = [];
            end
            
            obj.untrimmed_AA_centers = [xv,yv];
            scatter((obj.untrimmed_AA_centers(:,1)-1)*obj.scan_stepsize,(obj.untrimmed_AA_centers(:,2)-1)*obj.scan_stepsize,100,'b','filled');
            % ascertain proximity to array boundary
            [r,c] = size(fAA);
            adist = min([abs(1-yv),abs(r-yv),abs(1-xv),abs(c-xv)],[],2);
            xv(adist < edge_threshold) = [];
            yv(adist < edge_threshold) = [];
            AA_centroids = horzcat(xv,yv);
            scatter((AA_centroids(:,1)-1)*obj.scan_stepsize,(AA_centroids(:,2)-1)*obj.scan_stepsize,100,'r','filled');
        end
        
        
        
        
        
        
        function [meanAA_DP,fitting_values] = fitAAToCircle(obj,plot_average_AA,filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,use_robust_multistart,AA_amp_threshold,AAblur_Gaussian_sigma)
            if nargin < 2
                plot_average_AA = 0;
            end
            if nargin < 3
                filteramp = [];
                filterrange = [];
            end
            if nargin < 5
                hotspot_suppression = false;
            end
            if nargin < 6
                edge_threshold = 10;
            end
            if nargin < 7
                gauss_filter_threshold = 0.1;
            end
            if nargin < 8
                use_robust_multistart = 0;
            end
            if nargin < 9
                AA_amp_threshold = 0.71;
            end
            if nargin < 10
                AAblur_Gaussian_sigma = 5;
            end
            
            [AA_centroids,amask,f] = obj.findAAcenters(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma);
            nAA = size(AA_centroids,1);
            [r,c] = size(amask);
            xbase = 1:c;
            ybase = 1:r;
            [xspace,yspace] = meshgrid(xbase,ybase);
            options.Display = 'off';
            % Method 1: fit a circle to the thresholded region.
            predfun = @(x0,y0,r) isInCircle(xspace,yspace,x0,y0,r);
            fitstorage_circle = zeros(nAA,3);
            linind_all_storage = [];
            if use_robust_multistart == 1
                test_radii = 2:0.25:8; % Consider if we should introduce small perturbations to the com location or not.
            else
                test_radii = 5;
            end
            outer_error_storage = zeros(nAA,1);
            for i = 1:nAA
                fprintf('Fitting AA circle %d of %d.\n',i,nAA);
                thisAAcentroid = AA_centroids(i,:);
                thismask = amask;
                thismask(~isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),50)) = 0;
                fitfun = @(params) nnz(thismask - predfun(params(1),params(2),params(3)));
                % NPK added 04/25/2020: the inner multistart loop.
                error_storage = zeros(numel(test_radii),1);
                fitval_storage = zeros(numel(test_radii),3);
                if use_robust_multistart == 2 % Fit with PSwarm
                    center_increment = 4;  % location of center will probably not vary any more than that.
                    Problem.ObjFunction = fitfun;
                    Problem.LB = [thisAAcentroid-center_increment,0];
                    Problem.UB = [thisAAcentroid+center_increment,20];
                    Options.MaxObj = 2000;
                    Options.Size = 100;
                    Options.IPrint = -1;
                    [thisfitvals, error_storage] = PSwarm(Problem,[],Options);
                    thisfitvals = fminsearch(fitfun,thisfitvals,options);
                else % Fit with fminsearch
                    for j = 1:numel(test_radii)
                        init_guesses = [thisAAcentroid,test_radii(j)];
                        thisfitvals = fminsearch(fitfun,init_guesses,options);
                        error_storage(j) = fitfun(thisfitvals);
                        fitval_storage(j,:) = thisfitvals;
                    end
                end
                if use_robust_multistart == 1
                    [besterr,best_idx] = min(error_storage);
                    fitvals = fitval_storage(best_idx,:);
                    outer_error_storage(i) = besterr;
                else
                    fitvals = thisfitvals;
                    outer_error_storage(i) = error_storage;
                end
                fitted_circle = predfun(fitvals(1),fitvals(2),fitvals(3));
                linind_all_storage = [linind_all_storage;find(fitted_circle)];
                if false
                    figure;
                    imagesc(cat(3,thismask,zeros(size(thismask)),zeros(size(thismask))));
                    hold on
                    h = imagesc(cat(3,zeros(size(fitted_circle)),fitted_circle,zeros(size(fitted_circle))));
                    h.AlphaData = 0.5;
                    axis equal
                end
                fitstorage_circle(i,:) = fitvals;
            end % This does not appear to be terribly accurate -- it's OK
            fitstorage_circle(:,3) = abs(fitstorage_circle(:,3));  % because of the occasional negative radius problem
            fitting_values = fitstorage_circle;
            
            figure(f)
            hold on
            for i = 1:nAA
                cent = fitstorage_circle(i,1:2)*obj.scan_stepsize;
                %         cent2 = fitstorage_circle(i,2);
                radius = fitstorage_circle(i,3)*obj.scan_stepsize;
                h=viscircles(cent,radius,'Color','r');
            end
            axis equal
            
            % Get averaged AA diffraction patterns for Maddie
            if plot_average_AA
                runningDP = zeros([obj.datacube_size(1:2)]);
                npatterns = numel(linind_all_storage);
                for i = 1:npatterns
                    [to_use_subind1,to_use_subind2] = ind2sub(size(amask),linind_all_storage(i));
                    DP = obj.singleLoad(to_use_subind1,to_use_subind2);
                    runningDP = runningDP + double(DP);
                end
                meanAA_DP = runningDP./npatterns;
                figure;
                imagesc(meanAA_DP);
                colormap gray
                axis equal
                colorbar
                set(gca,'yDir','normal');
                obj.averageAA_DP = meanAA_DP;
            else
                meanAA_DP = [];
            end
            
            obj.AA_circlefit_values = fitstorage_circle;
            obj.AA_circlefit_pixelerrors = outer_error_storage;
        end
        
        
        
        
        
        
        % New function 06262020: bootstrap the circle fit. Requires
        % refactoring the fit so that it proceeds in linear fashion.
        % Arguments 2-5 should matched to the values used for the original
        % circlefit.
        function [boot_95CIs,boot_stderrs,bootstrap_bestcoords,bootstrap_bestresiduals] = getAAFitUncertainty(obj,gauss_filter_threshold,use_robust_multistart,AA_amp_threshold,AAblur_Gaussian_sigma,filteramp,filterrange,bidx,nboot)
            % (1) Choose AA region to bootstrap
            if nargin < 6
                stringcell = {'soft multistart','amplitude',[0.3,5]};
                filterstruct = buildFilterStruct(stringcell);
            end
            if nargin < 7
                bidx = [];
            end
            if isempty(bidx)
                obj.makeCustomDisplacementColorPlot([],[],filteramp,filterrange);
                % Need to convert xy pixels to xy nm for superimposition
                AA_nm = (obj.AA_circlefit_values-1).*obj.scan_stepsize;
                hold on
                scatter(AA_nm(:,1),AA_nm(:,2),20,'g','filled');
                disp('Select an AA region to bootstrap. Will snap to selection.');
                [coords] = ginput(1);
                xdiff = AA_nm(:,1) - coords(1);
                ydiff = AA_nm(:,2) - coords(2);
                dists = sqrt(xdiff.^2 + ydiff.^2);
                [~,bidx] = min(dists);
                scatter(AA_nm(bidx,1),AA_nm(bidx,2),50,'r','filled');
                fprintf('Index to bootstrap is %d.\n',bidx);
            end
            hotspot_suppression = 0;
            edge_threshold = 0;
            
            % Obtain thresholded data for the fit. This is accomplished by
            % amask.
            [AA_centroids,amask,f] = obj.findAAcenters(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma);
            nAA = size(AA_centroids,1);
            [r,c] = size(amask);
            % Now, crop a window and rework the fit into linear fashion for
            % bootstrapping. Will crop a 20x20 window around the original
            % AA fit.
            coords = obj.AA_circlefit_values(bidx,:);
            x0 = round(coords(1));
            y0 = round(coords(2));
            xbasecropped = (x0-10):(x0+10);
            ybasecropped = (y0-10):(y0+10);
            [xspacecropped,yspacecropped] = meshgrid(xbasecropped,ybasecropped);
            amask_cropped = amask(ybasecropped,xbasecropped);  
            % Coords are in xy pixels, so to access correct region through rc indexing, need to flip.
            
            % Verify we can reproduce fit.
            % Obtain linear indices for the data used for the fit, in
            % linear form.
            [cropcoords1,cropcoords2] = meshgrid(1:numel(xbasecropped),1:numel(ybasecropped));
            sz = size(cropcoords1);  
            spacecoords_lin = sub2ind(sz,cropcoords2(:),cropcoords1(:));
            fitdata_lin = amask_cropped(:);
            
            fitfun = @(c) AAmask_linearfitfun( c(1),c(2),c(3),xspacecropped,yspacecropped,spacecoords_lin,fitdata_lin );
            % initial optimization
%             linear_orgfit_coords = fminsearch(fitfun,round(coords));
            % Well, it didn't find the original best values, but it did
            % optimize to something similar.
            if use_robust_multistart == 1
                test_radii = 2:0.25:8; % Consider if we should introduce small perturbations to the com location or not.
                missvals = zeros(numel(test_radii),1);
                coordstor = zeros(numel(test_radii),3);
                for i = 1:numel(test_radii)
                    thisfit_coords = fminsearch(fitfun,[coords(1:2),test_radii(i)]);
                    missvals(i) = fitfun(thisfit_coords);
                    coordstor(i,:) = thisfit_coords;
                end
                [bestmissval,bestidx] = min(missvals);
                bestfitcoords = coordstor(bestidx,:);
            end
            % This fit is really lousy and unstable. But perhaps not all of
            % it is due to the uncertainty in the data -- some is just due
            % to the crappy nature of fitting a mask.
            
            
            % Perform the bootstrap
            bootmax = max(spacecoords_lin); % in linear indices
            bootstrap_bestcoords = zeros(nboot,3);
            bootstrap_bestresiduals = zeros(nboot,1);
            for b = 1:nboot
                if mod(b,10) == 0
                    fprintf('Bootstrap iteration %d of %d.\n',b,nboot);
                end
                % Randomly sample linear indices with replacement.
                bootinds = randi(bootmax,[bootmax,1]);  % because bootmax both gives the range of acceptible values and the total number to pull.
                bootdata = amask_cropped(bootinds);
                % xspacecropped and yspacecropped remain the same for
                % correct production of the prediction mask.
                % bootinds and bootdata are giving the linear indices of
                % that mask to draw from for the various pieces of data.
                boot_fitfun = @(c) AAmask_linearfitfun( c(1),c(2),c(3),xspacecropped,yspacecropped,bootinds,bootdata );
                % Perform multistart fitting
                if use_robust_multistart == 1
                    test_radii = 2:0.25:8; % Consider if we should introduce small perturbations to the com location or not.
                    missvals = zeros(numel(test_radii),1);
                    coordstor = zeros(numel(test_radii),3);
                    for i = 1:numel(test_radii)
                        thisfit_coords = fminsearch(boot_fitfun,[coords(1:2),test_radii(i)]);
                        missvals(i) = fitfun(thisfit_coords);
                        coordstor(i,:) = thisfit_coords;
                    end
                    [bestmissval,bestidx] = min(missvals);
                    bestfitcoords = coordstor(bestidx,:);
                    bootstrap_bestresiduals(b) = bestmissval;
                    bootstrap_bestcoords(b,:) = bestfitcoords;
                end
            end
            
            % Get the 95% confidence interval
            boot_95CIs = vertcat(prctile(bootstrap_bestcoords,2.5,1),prctile(bootstrap_bestcoords,97.5,1));
            % Get the standard error
            boot_stderrs = std(bootstrap_bestcoords);
        end
        
        
        
        
        
        % Requires having fit the AA regions to a circular mask via the
        % above function, fitAAToCircle
        function tf_boolmat = getAACircleMask(obj,radius_multiplier,min_radius,include_untrimmed)
            if isempty(obj.AA_circlefit_values)
                %                 obj.fitAAToCircle(1);
                pos = obj.AA_Gaussian_Circle_Fit_Params(:,1:2);
                rads = obj.AA_Gaussian_Circle_Fit_Params(:,3);
            else
                pos = obj.AA_circlefit_values(:,1:2);
                rads = obj.AA_circlefit_values(:,3);
            end
            if nargin < 2
                radius_multiplier = 1;
            end
            if (nargin < 3), min_radius = []; end;
            if (nargin < 4), include_untrimmed = false; end;
               
            
            if include_untrimmed
                pos = vertcat(pos,obj.untrimmed_AA_centers);
                rads = vertcat(rads,repmat(min_radius,size(obj.untrimmed_AA_centers,1),1));
            end
            
            xbase = 1:obj.datacube_size(4);
            ybase = 1:obj.datacube_size(3);
            [x,y] = meshgrid(xbase,ybase);
            tf_boolmat = false(size(x));
            radii_touse = rads;
            if ~isempty(min_radius)
                radii_touse(radii_touse < min_radius) = min_radius;
            end
            radii_touse = radii_touse*radius_multiplier;
            for i = 1:size(pos,1)
                tf = isInCircle(x,y,pos(i,1),pos(i,2),radii_touse(i));
                tf_boolmat = tf_boolmat | tf;
%                 tf_boolmat = [tf_boolmat;find(tf)];
            end
        end
        
        
        
        
        
        
        
        % New fitting function: Gaussian circle. There is some reason to
        % thing this might be the most stable of all the fits.
        function fitAAtoGaussianCircle(obj,filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,threshold_val,AAblur_Gaussian_sigma)
            % Basically use the exact same machinery as
            % fitAAtoGaussianEllipse() to obtain amplitude data for
            % fitting. Only change is which fitting functions are used, and
            % how the data save gets handled.
            if nargin < 8
                make_plots = false;
            end
            if nargin < 9
                mask_radius = 25;
            end
            if nargin < 11
                AAblur_Gaussian_sigma = 5;
            end
            if filteramp > 0
                displacement_field = obj.displacementContinuityFilter(filteramp,filterrange);
            else
                displacement_field = obj.DSC_fit_storage;
            end
            amplitude = sqrt(displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2);
            
            [AA_centroids,amask,f] = obj.findAAcenters(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma);
            nAA = size(AA_centroids,1);
            [r,c] = size(amask);
            xbase = 1:c;
            ybase = 1:r;
            [xspace,yspace] = meshgrid(xbase,ybase);
            % Lose semiminor and ellipse rotation parameters.
            fitstorage_Gcircle = zeros(nAA,6);
            % Set up the SP filter
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            [ ~, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
            H = HSV_color_stack(:,:,1);
            S = HSV_color_stack(:,:,2);
            V = HSV_color_stack(:,:,3);
            mask1 = S > 0.5;
            mask2 = V > 0.9;
            mask3 = mask1 & mask2;
            clear c
            
            % Begin loop to fit the AA domains
            for i = 1:nAA
                fprintf('Fitting circular Gaussian to AA domain %d of %d.\n',i,nAA);
                thisAAcentroid = AA_centroids(i,:);
                thismask = false(size(xspace));
                plotmask = false(size(xspace));
                thismask(isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),mask_radius)) = true;
                thismask(mask3) = false;
                xrectmaskcoords = thisAAcentroid(1)-mask_radius:thisAAcentroid(1)+mask_radius;
                yrectmaskcoords = thisAAcentroid(2)-mask_radius:thisAAcentroid(2)+mask_radius;
                xrectmaskcoords(xrectmaskcoords<1) = [];
                yrectmaskcoords(yrectmaskcoords<1) = [];
                xrectmaskcoords(xrectmaskcoords>size(xspace,1)) = [];
                yrectmaskcoords(yrectmaskcoords>size(xspace,2)) = [];
                plotmask(yrectmaskcoords,xrectmaskcoords) = true;
                %         residfun = @(c) ellipticGaussianResidfun( c, space, amplitude, thismask );%thismask(:).*(amplitude(:) - gausspredfun(c(1),c(2),c(3),c(4),c(5),c(6),c(7)));
                amplitude_to_fit = amplitude(thismask);
                xspace_touse = xspace(thismask);
                yspace_touse = yspace(thismask);
                space_tofit = horzcat(xspace_touse(:),yspace_touse(:));
                [fitvals,fitrmsr] = fitGaussianCircleInnards(obj,amplitude,amplitude_to_fit,space_tofit,thisAAcentroid,make_plots,plotmask,mask_radius,xrectmaskcoords,yrectmaskcoords,xspace,yspace);
                fitstorage_Gcircle(i,:) = [fitvals,fitrmsr];
            end
            
            % Nice thing is that we don't have to unpack the parameters
            % further. They are almost ready to go. Only need now to solve
            % the comparatively simple problem of figuring out what the
            % level curve equal to 0.71 Angstrom is, given the parameters
            % and the FWHM.
            final_AAparams = zeros(nAA,3);
            for i = 1:nAA
                % Convert FWHN (parameter 3), to radius at level curve.
                x0 = fitstorage_Gcircle(i,1);
                y0 = fitstorage_Gcircle(i,2);
                FWHM = fitstorage_Gcircle(i,3);
                A = fitstorage_Gcircle(i,4);
                B = fitstorage_Gcircle(i,5);
                sigma = FWHM/(2*sqrt(2*log(2)));
                r = sigma*sqrt(2*log(B/(A-0.71)));
                if imag(r) > 1e-10
                    final_AAparams(i,:) = [nan,nan,nan];
                else 
                    final_AAparams(i,:) = [x0,y0,r];
                end
            end
            
            % Store
            obj.AA_Gaussian_Circle_Fit_Params = final_AAparams;
        end
        
        
        
        
        
        % New fitting function: Gaussian circle. There is some reason to
        % thing this might be the most stable of all the fits.
        function [stderrs_of_fit,CI_lb,CI_ub,raw_bootstrap_data] = bootstrapAAGaussianCircle(obj,filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,threshold_val,AAblur_Gaussian_sigma,bootstrap_idx,nboot)
            % Basically use the exact same machinery as
            % fitAAtoGaussianEllipse() to obtain amplitude data for
            % fitting. Only change is which fitting functions are used, and
            % how the data save gets handled.
            if nargin < 8
                make_plots = false;
            end
            if nargin < 9
                mask_radius = 25;
            end
            if nargin < 11
                AAblur_Gaussian_sigma = 5;
            end
            if filteramp > 0
                displacement_field = obj.displacementContinuityFilter(filteramp,filterrange);
            else
                displacement_field = obj.DSC_fit_storage;
            end
            amplitude = sqrt(displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2);
            
            [AA_centroids,amask,f] = obj.findAAcenters(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma);
            nAA = size(AA_centroids,1);
            [r,c] = size(amask);
            xbase = 1:c;
            ybase = 1:r;
            [xspace,yspace] = meshgrid(xbase,ybase);
            % Lose semiminor and ellipse rotation parameters.
            fitstorage_Gcircle = zeros(nAA,6);
            % Set up the SP filter
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            [ ~, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
            H = HSV_color_stack(:,:,1);
            S = HSV_color_stack(:,:,2);
            V = HSV_color_stack(:,:,3);
            mask1 = S > 0.5;
            mask2 = V > 0.9;
            mask3 = mask1 & mask2;

            clear c
            % Changed up to get the bootstrapping going.
            % The only two variables that actually go into the fitting function are space_tofit
            % and amplitude_to_fit. Therefore bootstrap these, but do not bootstrap all this other stuff.
            thisAAcentroid = AA_centroids(bootstrap_idx,:);
            thismask = false(size(xspace));
            plotmask = false(size(xspace));
            thismask(isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),mask_radius)) = true;
            thismask(mask3) = false;
            xrectmaskcoords = thisAAcentroid(1)-mask_radius:thisAAcentroid(1)+mask_radius;
            yrectmaskcoords = thisAAcentroid(2)-mask_radius:thisAAcentroid(2)+mask_radius;
            xrectmaskcoords(xrectmaskcoords<1) = [];
            yrectmaskcoords(yrectmaskcoords<1) = [];
            xrectmaskcoords(xrectmaskcoords>size(xspace,1)) = [];
            yrectmaskcoords(yrectmaskcoords>size(xspace,2)) = [];
            plotmask(yrectmaskcoords,xrectmaskcoords) = true;
            %         residfun = @(c) ellipticGaussianResidfun( c, space, amplitude, thismask );%thismask(:).*(amplitude(:) - gausspredfun(c(1),c(2),c(3),c(4),c(5),c(6),c(7)));
            amplitude_to_fit_org = amplitude(thismask);
            xspace_touse = xspace(thismask);
            yspace_touse = yspace(thismask);
            space_tofit_org = horzcat(xspace_touse(:),yspace_touse(:));
            % bootstrap rows
            % Bootstrap
            nrows = size(amplitude_to_fit_org,1);
            boot_fitstorage = zeros(nboot,6);
            for i = 1:nboot
                if mod(i,10) == 0
                    fprintf('AA Gaussian ellipse bootstrap iteration %d of %d.\n',i,nboot);
                end
                bootinds = randi(nrows,[nrows,1]);
                amplitude_to_fit_boot = amplitude_to_fit_org(bootinds,:);
                space_tofit_boot = space_tofit_org(bootinds,:);
                [fitvals,fitrmsr] = fitGaussianCircleInnards(obj,amplitude,amplitude_to_fit_boot,space_tofit_boot,thisAAcentroid,make_plots,plotmask,mask_radius,xrectmaskcoords,yrectmaskcoords,xspace,yspace);
                boot_fitstorage(i,:) = [fitvals,fitrmsr];
            end
            fitstorage_Gcircle = boot_fitstorage;
            
            % Nice thing is that we don't have to unpack the parameters
            % further. They are almost ready to go. Only need now to solve
            % the comparatively simple problem of figuring out what the
            % level curve equal to 0.71 Angstrom is, given the parameters
            % and the FWHM.
            final_AAparams = zeros(nAA,3);
            for i = 1:nboot
                % Convert FWHN (parameter 3), to radius at level curve.
                x0 = fitstorage_Gcircle(i,1);
                y0 = fitstorage_Gcircle(i,2);
                FWHM = fitstorage_Gcircle(i,3);
                A = fitstorage_Gcircle(i,4);
                B = fitstorage_Gcircle(i,5);
                sigma = FWHM/(2*sqrt(2*log(2)));
                r = sigma*sqrt(2*log(B/(A-0.71)));
                if imag(r) > 1e-10
                    final_AAparams(i,:) = [nan,nan,nan];
                else 
                    final_AAparams(i,:) = [x0,y0,r];
                end
            end
            
            % Compute standard errors, percentile confidence intervals, and
            % return.
            stderrs_of_fit = std(final_AAparams);
            CI_lb = prctile(final_AAparams,2.5,1);
            CI_ub = prctile(final_AAparams,97.5,1);
            raw_bootstrap_data = final_AAparams;
        end
        
        
        
        
        
        
        function [radius_mean_nm,radius_std_nm,radius_stderr,radius_95confdelta] = getGaussianCircleFitStatistics(obj,radiistatsrange)
            params_to_use = obj.AA_Gaussian_Circle_Fit_Params;
            radii = params_to_use(:,3)*obj.scan_stepsize;
            if nargin >= 2
                remove_idxs = min(radiistatsrange) > radii | max(radiistatsrange) < radii;
                radii(remove_idxs,:) = [];
                params_to_use(remove_idxs,:) = [];
            end
            try 
                obj.makeCustomDisplacementColorPlot([],[],0.3,5,1,0,0,0);
            catch
                % Multistart filter not defined for simulated data.
                obj.makeCustomDisplacementColorPlot([],[],[],[],1,0,0,0);
            end
            hold on
            centsplot = (params_to_use(:,1:2) - 1)*obj.scan_stepsize;
            
            scatter(centsplot(:,1),centsplot(:,2),'g','filled');
            viscircles(centsplot,radii,'Color','g');
            
            
            
            radius_mean_nm = mean(radii,'omitnan');
            radius_std_nm = std(radii,'omitnan');
            N = size(params_to_use,1);
            radius_stderr = radius_std_nm/sqrt(N);
            radius_95confdelta = 2*radius_stderr;
            figure;
            histogram(radii);
            title('AA Gaussian circle fit radii');
            xlabel('Radius (nm)');
            ylabel('Counts');
        end
        
        
        
        
        
        
        % Copied over from fitGaussianEllipseInnards. Just change the
        % fitting function definitions on the inside. This sort of fit
        % should be quite robust -- don't need multistart.
        function [fitvals,fitrmsr] = fitGaussianCircleInnards(obj,amplitude,amplitude_to_fit,space_tofit,thisAAcentroid,make_plots,plotmask,mask_radius,xrectmaskcoords,yrectmaskcoords,xspace,yspace)
            options = optimoptions('lsqnonlin');
            options.MaxFunEvals = 20000;
            options.MaxIter = 2000;
            options.CheckGradients = false;
            options.SpecifyObjectiveGradient = false;
            options.Display = 'off';
            
            
%             gausspredfun = @(c) ellipticGaussianPredfun2( space_tofit,c );
            residfun = @(c) circularGaussianResidfun2( c,space_tofit,amplitude_to_fit );
            rmsrfun = @(c) rms(residfun(c));
            
            %         thisAAcentroid = fliplr(thisAAcentroid);
            init_guess = [thisAAcentroid,5,1.4,1];
            lb = [thisAAcentroid-mask_radius,0,0,0];
            ub = [thisAAcentroid+mask_radius,100,5,5];
            fitvals = lsqnonlin(residfun,init_guess,lb,ub,options);
            fitrmsr = rmsrfun(fitvals);
            
            amplitude_to_plot = amplitude(plotmask);
            amplitude_to_plug_in = amplitude(xrectmaskcoords,yrectmaskcoords);
            xspace_toplot = xspace(plotmask);
            yspace_toplot = yspace(plotmask);
            space_toplot = horzcat(xspace_toplot(:),yspace_toplot(:));
            fitted_gaussian_toplot = circularGaussianPredfun2( space_toplot,fitvals );
            fitted_gaussian_toplot = reshape(fitted_gaussian_toplot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
            amplitude_to_plot = reshape(amplitude_to_plot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
            if make_plots
                thismask_cropped = thismask(yrectmaskcoords,xrectmaskcoords);
                amplitude_to_plot2 = amplitude_to_plot;
                amplitude_to_plot2(~thismask_cropped) = nan;
                figure
                pcolor(xrectmaskcoords,yrectmaskcoords,-amplitude_to_plot2);
                axis equal
                shading flat
                title 'Amplitude data to fit'
                xlabel('xcoords');
                ylabel('ycoords');
                colorbar
                
                fitted_gaussian_toplot(~thismask_cropped) = nan;
                figure
                pcolor(xrectmaskcoords,yrectmaskcoords,-fitted_gaussian_toplot);
                axis equal
                shading flat
                title 'Gaussian fit'
                xlabel('xcoords');
                ylabel('ycoords');
                colorbar
                
                resid_toplot = circularGaussianResidfun2( fitvals,space_toplot,amplitude_to_plot(:) );
                resid_toplot = reshape(resid_toplot,size(thismask_cropped));
                resid_toplot(~thismask_cropped) = nan;
                figure
                pcolor(xrectmaskcoords,yrectmaskcoords,resid_toplot);
                axis equal
                shading flat
                title 'Residuals'
                xlabel('xcoords');
                ylabel('ycoords');
                colorbar
            end
        end
        
        
        
        
        
        
        
        % Refactored out of fit_AA_to_ellipse
        % Unlike the fitAAToCircle method, this function only uses the
        % gauss_filter_threshold variable to get a first approximation for
        % the locations of the AA regions. This will be used to start the
        % Gaussian fit, which is a non-thresholded process.
        function converge_storage = fitAAtoGaussianEllipse(obj,filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,threshold_val,AAblur_Gaussian_sigma)
            if nargin < 8
                make_plots = false;
            end
            if nargin < 9
                mask_radius = 25;
            end
            if nargin < 11
                AAblur_Gaussian_sigma = 5;
            end
            if filteramp > 0
                displacement_field = obj.displacementContinuityFilter(filteramp,filterrange);
            else
                displacement_field = obj.DSC_fit_storage;
            end
            amplitude = sqrt(displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2);
            
            [AA_centroids,amask,f] = obj.findAAcenters(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma);
            nAA = size(AA_centroids,1);
            [r,c] = size(amask);
            xbase = 1:c;
            ybase = 1:r;
            [xspace,yspace] = meshgrid(xbase,ybase);
            fitstorage_2DGauss = zeros(nAA,8);
            % Set up the SP filter
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            [ ~, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
            H = HSV_color_stack(:,:,1);
            S = HSV_color_stack(:,:,2);
            V = HSV_color_stack(:,:,3);
            mask1 = S > 0.5;
            mask2 = V > 0.9;
            mask3 = mask1 & mask2;
            
            
            
            clear c
            converge_storage = zeros(nAA,1);
            for i = 1:nAA
                fprintf('Fitting elliptical Gaussian to AA domain %d of %d.\n',i,nAA);
                thisAAcentroid = AA_centroids(i,:);
                thismask = false(size(xspace));
                plotmask = false(size(xspace));
                thismask(isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),mask_radius)) = true;
                thismask(mask3) = false;
                xrectmaskcoords = thisAAcentroid(1)-mask_radius:thisAAcentroid(1)+mask_radius;
                yrectmaskcoords = thisAAcentroid(2)-mask_radius:thisAAcentroid(2)+mask_radius;
                xrectmaskcoords(xrectmaskcoords<1) = [];
                yrectmaskcoords(yrectmaskcoords<1) = [];
                xrectmaskcoords(xrectmaskcoords>size(xspace,1)) = [];
                yrectmaskcoords(yrectmaskcoords>size(xspace,2)) = [];
                plotmask(yrectmaskcoords,xrectmaskcoords) = true;
                %         residfun = @(c) ellipticGaussianResidfun( c, space, amplitude, thismask );%thismask(:).*(amplitude(:) - gausspredfun(c(1),c(2),c(3),c(4),c(5),c(6),c(7)));
                amplitude_to_fit = amplitude(thismask);
                xspace_touse = xspace(thismask);
                yspace_touse = yspace(thismask);
                space_tofit = horzcat(xspace_touse(:),yspace_touse(:));
                [fitvals,fitrmsr,num_global_converged] = fitGaussianEllipseInnards(obj,amplitude,amplitude_to_fit,space_tofit,thisAAcentroid,make_plots,plotmask,mask_radius,xrectmaskcoords,yrectmaskcoords,xspace,yspace);
                converge_storage(i) = num_global_converged;
                % beginning of refactored region
                %                 thisAAcentroid = AA_centroids(i,:);
                %                 thismask = false(size(xspace));
                %                 plotmask = false(size(xspace));
                %                 thismask(isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),mask_radius)) = true;
%                 thismask(mask3) = false;
%                 xrectmaskcoords = thisAAcentroid(1)-mask_radius:thisAAcentroid(1)+mask_radius;
%                 yrectmaskcoords = thisAAcentroid(2)-mask_radius:thisAAcentroid(2)+mask_radius;
%                 xrectmaskcoords(xrectmaskcoords<1) = [];
%                 yrectmaskcoords(yrectmaskcoords<1) = [];
%                 xrectmaskcoords(xrectmaskcoords>size(xspace,1)) = [];
%                 yrectmaskcoords(yrectmaskcoords>size(xspace,2)) = [];
%                 plotmask(yrectmaskcoords,xrectmaskcoords) = true;
%                 %         residfun = @(c) ellipticGaussianResidfun( c, space, amplitude, thismask );%thismask(:).*(amplitude(:) - gausspredfun(c(1),c(2),c(3),c(4),c(5),c(6),c(7)));
%                 amplitude_to_fit = amplitude(thismask);
%                 xspace_touse = xspace(thismask);
%                 yspace_touse = yspace(thismask);
%                 space_tofit = horzcat(xspace_touse(:),yspace_touse(:));
%                 gausspredfun = @(c) ellipticGaussianPredfun2( space_tofit,c );
%                 residfun = @(c) ellipticGaussianResidfun2( c,space_tofit,amplitude_to_fit );
%                 rmsrfun = @(c) rms(residfun(c));
%                 
%                 %         thisAAcentroid = fliplr(thisAAcentroid);
%                 init_guesses = [thisAAcentroid,4,4,0,1.4,1;
%                                 thisAAcentroid,10,2,0,1.2,1;
%                                 thisAAcentroid,2,10,0,1.2,1;
%                                 thisAAcentroid,6,6,0.8,1.2,1;
%                                 thisAAcentroid,6,6,-0.8,1.2,1;
%                                 thisAAcentroid+[2,0],6,6,0,1.2,1;
%                                 thisAAcentroid+[0,2],6,6,0,1.2,1;
%                                 thisAAcentroid-[2,0],6,6,0,1.2,1;
%                                 thisAAcentroid-[0,2],6,6,0,1.2,1];
%                 lb = [thisAAcentroid-mask_radius,0,0,-1,0,0];
%                 ub = [thisAAcentroid+mask_radius,100,100,1,5,5];
%                 
%                 fitvals_stor = zeros(size(init_guesses));
%                 fitrmsr_stor = zeros(size(init_guesses,1),1);
%                 for q = 1:size(init_guesses,1)
%                     fitvals_stor(q,:) = lsqnonlin(residfun,init_guesses(q,:),lb,ub,options);
%                     fitrmsr_stor(q) = rmsrfun(fitvals_stor(q,:));
%                 end
%                 
%                 roundconstant = 5;
%                 [bestval,bestidx] = min(fitrmsr_stor);
%                 roundrmsrstor = round(fitrmsr_stor,roundconstant);
%                 num_global_converged = nnz(roundrmsrstor == round(bestval,roundconstant));
%                 fitvals = fitvals_stor(bestidx,:);
%                 fitrmsr = fitrmsr_stor(bestidx);
%                 converge_storage(i) = num_global_converged;
%                 
%                 amplitude_to_plot = amplitude(plotmask);
%                 amplitude_to_plug_in = amplitude(xrectmaskcoords,yrectmaskcoords);
%                 xspace_toplot = xspace(plotmask);
%                 yspace_toplot = yspace(plotmask);
%                 space_toplot = horzcat(xspace_toplot(:),yspace_toplot(:));
%                 fitted_gaussian_toplot = ellipticGaussianPredfun2( space_toplot,fitvals );
%                 %         fitted_gaussian = gausspredfun(fitvals);
%                 fitted_gaussian_toplot = reshape(fitted_gaussian_toplot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
%                 amplitude_to_plot = reshape(amplitude_to_plot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
%                 if make_plots
%                     %             figure;
%                     %             imagesc(cat(3,amplitude,zeros(size(amplitude_to_plot)),zeros(size(amplitude_to_plot))));
%                     %             hold on
%                     %             h = imagesc(cat(3,zeros(size(fitted_gaussian)),fitted_gaussian,zeros(size(fitted_gaussian))));
%                     %             h.AlphaData = 0.5;
%                     %             axis equal
%                     thismask_cropped = thismask(yrectmaskcoords,xrectmaskcoords);
%                     amplitude_to_plot2 = amplitude_to_plot;
%                     amplitude_to_plot2(~thismask_cropped) = nan;
%                     figure
%                     pcolor(xrectmaskcoords,yrectmaskcoords,-amplitude_to_plot2);
%                     axis equal
%                     shading flat
%                     title 'Amplitude data to fit'
%                     xlabel('xcoords');
%                     ylabel('ycoords');
%                     colorbar
%                     
%                     fitted_gaussian_toplot(~thismask_cropped) = nan;
%                     figure
%                     pcolor(xrectmaskcoords,yrectmaskcoords,-fitted_gaussian_toplot);
%                     axis equal
%                     shading flat
%                     title 'Gaussian fit'
%                     xlabel('xcoords');
%                     ylabel('ycoords');
%                     colorbar
%                     
%                     resid_toplot = ellipticGaussianResidfun2( fitvals,space_toplot,amplitude_to_plot(:) );
%                     resid_toplot = reshape(resid_toplot,size(thismask_cropped));
%                     resid_toplot(~thismask_cropped) = nan;
%                     figure
%                     pcolor(xrectmaskcoords,yrectmaskcoords,resid_toplot);
%                     axis equal
%                     shading flat
%                     title 'Residuals'
%                     xlabel('xcoords');
%                     ylabel('ycoords');
%                     colorbar
%                 end
                % end of refactored region
                fitstorage_2DGauss(i,:) = [fitvals,fitrmsr];
            end
            fitstorage_2DGauss
            
            if make_plots
                figure(f)
                hold on
            end
            % Regardless of whether we are making the figure, we do still
            % need to use the code to extract FWHM and angle for storage.
            angle_axis_storage = zeros(nAA,5);  % let's put x0, y0 in as the fourth and fifth values
            for i = 1:nAA
                fprintf('Extracting FWHM parameters for AA domain %d of %d.\n',i,nAA);
                % fsolve to back out semimajor/minor axes and rotation angle
                FWHM1 = fitstorage_2DGauss(i,3);
                FWHM2 = fitstorage_2DGauss(i,4);
                p = fitstorage_2DGauss(i,5);
                s1 = FWHM1/(2*sqrt(2*log(2)));
                s2 = FWHM2/(2*sqrt(2*log(2)));
                if threshold_val == 'a'  % which stands for "all"
                    k = 2*log(2)*(1-p^2);
                else
                    A = fitstorage_2DGauss(i,6);
                    B = fitstorage_2DGauss(i,7);
                    k = 2*(1-p^2)*log(B/(A-threshold_val));
                end
                eqn1 = @(alpha,a,b) cos(alpha)^2/(a^2) + sin(alpha)^2/(b^2) - 1/(k*s1^2);
                eqn2 = @(alpha,a,b) sin(2*alpha)*(1/(a^2) - 1/(b^2)) + 2*p/(k*s1*s2);
                eqn3 = @(alpha,a,b) sin(alpha)^2/(a^2) + cos(alpha)^2/(b^2) - 1/(k*s2^2);
                eqns = @(x) [eqn1(x(1),x(2),x(3)),eqn2(x(1),x(2),x(3)),eqn3(x(1),x(2),x(3))];
                x0 = [p,FWHM1,FWHM2];
                options = optimset('lsqnonlin');
                options.TolFun = 1e-10;
                options.Display = 'off';
                % set bounds because there seems to be some sort of ambiguity left
                % in these equations
                if p > 0
                    lb = [0,0,0];
                    ub = [pi/2,inf,inf];
                else
                    lb = [-pi/2,0,0];
                    ub = [0,inf,inf];
                end
                solvevals = lsqnonlin(eqns,x0,lb,ub,options);
                residuals = eqns(solvevals);
                whilecount = 1;
                newx0s = [2*p,FWHM2,FWHM1;
                          lb(1),FWHM1,FWHM2;
                          ub(1),FWHM1,FWHM2;
                          p,2,10;
                          p,10,2];  % shake things up a little in the initial guess
                while rms(residuals) > 1e-5  % cycle through other initial guess possibilities if the equations are not solved well.
                    disp('Retrying fit');
                    newx0 = newx0s(whilecount,:);
%                     options.Display = 'final';
                    solvevals = lsqnonlin(eqns,newx0,lb,ub,options);
                    residuals = eqns(solvevals);
                    whilecount = whilecount + 1;
                    if whilecount > size(newx0s,1)  % overflow about to happen
                        disp('Parameter extraction not within tolerances in fitAAtoGaussianEllipse(). Check the quality of the fit.');
                        warning('Parameter extraction not within tolerances in fitAAtoGaussianEllipse(). Check the quality of the fit.');
                        break
                    end
                end
                
                semimajor = max(solvevals(2:3));
                semiminor = min(solvevals(2:3));
                angle = solvevals(1);  % no attempt here to phase wrap, as we don't really care about the angle value anyways.
                
                %         ang = fitstorage_2DGauss(i,7);
                %         ra = fitstorage_2DGauss(i,3);
                %         rb = fitstorage_2DGauss(i,4);
                % these are correct
                x0 = fitstorage_2DGauss(i,1);
                y0 = fitstorage_2DGauss(i,2);
                angle_axis_storage(i,:) = [angle,semimajor,semiminor,x0,y0];
                if make_plots
                    h=ellipse(semimajor*obj.scan_stepsize,semiminor*obj.scan_stepsize,angle,x0*obj.scan_stepsize,y0*obj.scan_stepsize);
                    axis equal
                end
            end
            obj.AAgaussianfit_FWHM_ellipses = angle_axis_storage;
            obj.AAgaussianfit_raw_parameters = fitstorage_2DGauss;
        end
        
        
        
        
        
             
        
        function [stderrs_of_fit,CI_lb,CI_ub,raw_bootstrap_data] = bootstrapAAGaussianEllipse(obj,filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,make_plots,mask_radius,threshold_val,AAblur_Gaussian_sigma,bootstrap_idx,nboot)
            if nargin < 8
                make_plots = false;
            end
            if nargin < 9
                mask_radius = 25;
            end
            if nargin < 11
                AAblur_Gaussian_sigma = 5;
            end
            if filteramp > 0
                displacement_field = obj.displacementContinuityFilter(filteramp,filterrange);
            else
                displacement_field = obj.DSC_fit_storage;
            end
            amplitude = sqrt(displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2);
            
            [AA_centroids,amask,f] = obj.findAAcenters(filteramp,filterrange,hotspot_suppression,edge_threshold,gauss_filter_threshold,AA_amp_threshold,AAblur_Gaussian_sigma);
            nAA = size(AA_centroids,1);
            [r,c] = size(amask);
            xbase = 1:c;
            ybase = 1:r;
            [xspace,yspace] = meshgrid(xbase,ybase);
            fitstorage_2DGauss = zeros(nAA,8);
            % Set up the SP filter
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            [ ~, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
            H = HSV_color_stack(:,:,1);
            S = HSV_color_stack(:,:,2);
            V = HSV_color_stack(:,:,3);
            mask1 = S > 0.5;
            mask2 = V > 0.9;
            mask3 = mask1 & mask2;
            
            
            
            clear c
            converge_storage = zeros(nAA,1);
            % Changed up to get the bootstrapping going.
            % The only two variables that actually go into the fitting function are space_tofit
            % and amplitude_to_fit. Therefore bootstrap these, but do not bootstrap all this other stuff.
            thisAAcentroid = AA_centroids(bootstrap_idx,:);
            thismask = false(size(xspace));
            plotmask = false(size(xspace));
            thismask(isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),mask_radius)) = true;
            thismask(mask3) = false;
            xrectmaskcoords = thisAAcentroid(1)-mask_radius:thisAAcentroid(1)+mask_radius;
            yrectmaskcoords = thisAAcentroid(2)-mask_radius:thisAAcentroid(2)+mask_radius;
            xrectmaskcoords(xrectmaskcoords<1) = [];
            yrectmaskcoords(yrectmaskcoords<1) = [];
            xrectmaskcoords(xrectmaskcoords>size(xspace,1)) = [];
            yrectmaskcoords(yrectmaskcoords>size(xspace,2)) = [];
            plotmask(yrectmaskcoords,xrectmaskcoords) = true;
            %         residfun = @(c) ellipticGaussianResidfun( c, space, amplitude, thismask );%thismask(:).*(amplitude(:) - gausspredfun(c(1),c(2),c(3),c(4),c(5),c(6),c(7)));
            amplitude_to_fit_org = amplitude(thismask);
            xspace_touse = xspace(thismask);
            yspace_touse = yspace(thismask);
            space_tofit_org = horzcat(xspace_touse(:),yspace_touse(:));
            % bootstrap rows
            nrows = size(amplitude_to_fit_org,1);
            boot_fitstorage = zeros(nboot,8);
            for i = 1:nboot
                if mod(i,10) == 0
                    fprintf('AA Gaussian ellipse bootstrap iteration %d of %d.\n',i,nboot);
                end
                bootinds = randi(nrows,[nrows,1]);
                amplitude_to_fit_boot = amplitude_to_fit_org(bootinds,:);
                space_tofit_boot = space_tofit_org(bootinds,:);
                [fitvals,fitrmsr,num_global_converged] = fitGaussianEllipseInnards(obj,amplitude,amplitude_to_fit_boot,space_tofit_boot,thisAAcentroid,make_plots,plotmask,mask_radius,xrectmaskcoords,yrectmaskcoords,xspace,yspace);
                converge_storage(i) = num_global_converged;
                
                boot_fitstorage(i,:) = [fitvals,fitrmsr];
            end
            
            if make_plots
                figure(f)
                hold on
            end
            % Regardless of whether we are making the figure, we do still
            % need to use the code to extract FWHM and angle for storage.
            angle_axis_storage = zeros(nboot,5);  % let's put x0, y0 in as the fourth and fifth values
            for i = 1:nboot
                fprintf('Extracting FWHM parameters for bootstrap iteration %d of %d.\n',i,nAA);
                % fsolve to back out semimajor/minor axes and rotation angle
                FWHM1 = boot_fitstorage(i,3);
                FWHM2 = boot_fitstorage(i,4);
                p = boot_fitstorage(i,5);
                s1 = FWHM1/(2*sqrt(2*log(2)));
                s2 = FWHM2/(2*sqrt(2*log(2)));
                if threshold_val == 'a'  % which stands for "all"
                    k = 2*log(2)*(1-p^2);
                else
                    A = boot_fitstorage(i,6);
                    B = boot_fitstorage(i,7);
                    k = 2*(1-p^2)*log(B/(A-threshold_val));
                end
                eqn1 = @(alpha,a,b) cos(alpha)^2/(a^2) + sin(alpha)^2/(b^2) - 1/(k*s1^2);
                eqn2 = @(alpha,a,b) sin(2*alpha)*(1/(a^2) - 1/(b^2)) + 2*p/(k*s1*s2);
                eqn3 = @(alpha,a,b) sin(alpha)^2/(a^2) + cos(alpha)^2/(b^2) - 1/(k*s2^2);
                eqns = @(x) [eqn1(x(1),x(2),x(3)),eqn2(x(1),x(2),x(3)),eqn3(x(1),x(2),x(3))];
                x0 = [p,FWHM1,FWHM2];
                options = optimset('lsqnonlin');
                options.TolFun = 1e-10;
                options.Display = 'off';
                % set bounds because there seems to be some sort of ambiguity left
                % in these equations
                if p > 0
                    lb = [0,0,0];
                    ub = [pi/2,inf,inf];
                else
                    lb = [-pi/2,0,0];
                    ub = [0,inf,inf];
                end
                solvevals = lsqnonlin(eqns,x0,lb,ub,options);
                residuals = eqns(solvevals);
                whilecount = 1;
                newx0s = [2*p,FWHM2,FWHM1;
                          lb(1),FWHM1,FWHM2;
                          ub(1),FWHM1,FWHM2;
                          p,2,10;
                          p,10,2];  % shake things up a little in the initial guess
                while rms(residuals) > 1e-5  % cycle through other initial guess possibilities if the equations are not solved well.
                    disp('Retrying fit');
                    newx0 = newx0s(whilecount,:);
%                     options.Display = 'final';
                    solvevals = lsqnonlin(eqns,newx0,lb,ub,options);
                    residuals = eqns(solvevals);
                    whilecount = whilecount + 1;
                    if whilecount > size(newx0s,1)  % overflow about to happen
                        disp('Parameter extraction not within tolerances in fitAAtoGaussianEllipse(). Check the quality of the fit.');
                        warning('Parameter extraction not within tolerances in fitAAtoGaussianEllipse(). Check the quality of the fit.');
                        break
                    end
                end
                
                semimajor = max(solvevals(2:3));
                semiminor = min(solvevals(2:3));
                angle = solvevals(1);  % no attempt here to phase wrap, as we don't really care about the angle value anyways.
                
                %         ang = fitstorage_2DGauss(i,7);
                %         ra = fitstorage_2DGauss(i,3);
                %         rb = fitstorage_2DGauss(i,4);
                % these are correct
                x0 = boot_fitstorage(i,1);
                y0 = boot_fitstorage(i,2);
                angle_axis_storage(i,:) = [angle,semimajor,semiminor,x0,y0];
                if make_plots
                    h=ellipse(semimajor*obj.scan_stepsize,semiminor*obj.scan_stepsize,angle,x0*obj.scan_stepsize,y0*obj.scan_stepsize);
                    axis equal
                end
            end
            
            % Compute standard errors, percentile confidence intervals, and
            % return.
            stderrs_of_fit = std(angle_axis_storage);
            CI_lb = prctile(angle_axis_storage,2.5,1);
            CI_ub = prctile(angle_axis_storage,97.5,1);
            raw_bootstrap_data = angle_axis_storage;
        end
        
        
        
        
        
        function [fitvals,fitrmsr,num_global_converged] = fitGaussianEllipseInnards(obj,amplitude,amplitude_to_fit,space_tofit,thisAAcentroid,make_plots,plotmask,mask_radius,xrectmaskcoords,yrectmaskcoords,xspace,yspace)
            options = optimoptions('lsqnonlin');
            options.MaxFunEvals = 20000;
            options.MaxIter = 2000;
            options.CheckGradients = false;
            options.SpecifyObjectiveGradient = false;
            options.Display = 'off';
            
            
            gausspredfun = @(c) ellipticGaussianPredfun2( space_tofit,c );
            residfun = @(c) ellipticGaussianResidfun2( c,space_tofit,amplitude_to_fit );
            rmsrfun = @(c) rms(residfun(c));
            
            %         thisAAcentroid = fliplr(thisAAcentroid);
            init_guesses = [thisAAcentroid,4,4,0,1.4,1;
                thisAAcentroid,10,2,0,1.2,1;
                thisAAcentroid,2,10,0,1.2,1;
                thisAAcentroid,6,6,0.8,1.2,1;
                thisAAcentroid,6,6,-0.8,1.2,1;
                thisAAcentroid+[2,0],6,6,0,1.2,1;
                thisAAcentroid+[0,2],6,6,0,1.2,1;
                thisAAcentroid-[2,0],6,6,0,1.2,1;
                thisAAcentroid-[0,2],6,6,0,1.2,1];
            lb = [thisAAcentroid-mask_radius,0,0,-1,0,0];
            ub = [thisAAcentroid+mask_radius,100,100,1,5,5];
            
            fitvals_stor = zeros(size(init_guesses));
            fitrmsr_stor = zeros(size(init_guesses,1),1);
            for q = 1:size(init_guesses,1)
                fitvals_stor(q,:) = lsqnonlin(residfun,init_guesses(q,:),lb,ub,options);
                fitrmsr_stor(q) = rmsrfun(fitvals_stor(q,:));
            end
            
            roundconstant = 5;
            [bestval,bestidx] = min(fitrmsr_stor);
            roundrmsrstor = round(fitrmsr_stor,roundconstant);
            num_global_converged = nnz(roundrmsrstor == round(bestval,roundconstant));
            fitvals = fitvals_stor(bestidx,:);
            fitrmsr = fitrmsr_stor(bestidx);
            
            amplitude_to_plot = amplitude(plotmask);
            amplitude_to_plug_in = amplitude(xrectmaskcoords,yrectmaskcoords);
            xspace_toplot = xspace(plotmask);
            yspace_toplot = yspace(plotmask);
            space_toplot = horzcat(xspace_toplot(:),yspace_toplot(:));
            fitted_gaussian_toplot = ellipticGaussianPredfun2( space_toplot,fitvals );
            %         fitted_gaussian = gausspredfun(fitvals);
            fitted_gaussian_toplot = reshape(fitted_gaussian_toplot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
            amplitude_to_plot = reshape(amplitude_to_plot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
            if make_plots
                %             figure;
                %             imagesc(cat(3,amplitude,zeros(size(amplitude_to_plot)),zeros(size(amplitude_to_plot))));
                %             hold on
                %             h = imagesc(cat(3,zeros(size(fitted_gaussian)),fitted_gaussian,zeros(size(fitted_gaussian))));
                %             h.AlphaData = 0.5;
                %             axis equal
                thismask_cropped = thismask(yrectmaskcoords,xrectmaskcoords);
                amplitude_to_plot2 = amplitude_to_plot;
                amplitude_to_plot2(~thismask_cropped) = nan;
                figure
                pcolor(xrectmaskcoords,yrectmaskcoords,-amplitude_to_plot2);
                axis equal
                shading flat
                title 'Amplitude data to fit'
                xlabel('xcoords');
                ylabel('ycoords');
                colorbar
                
                fitted_gaussian_toplot(~thismask_cropped) = nan;
                figure
                pcolor(xrectmaskcoords,yrectmaskcoords,-fitted_gaussian_toplot);
                axis equal
                shading flat
                title 'Gaussian fit'
                xlabel('xcoords');
                ylabel('ycoords');
                colorbar
                
                resid_toplot = ellipticGaussianResidfun2( fitvals,space_toplot,amplitude_to_plot(:) );
                resid_toplot = reshape(resid_toplot,size(thismask_cropped));
                resid_toplot(~thismask_cropped) = nan;
                figure
                pcolor(xrectmaskcoords,yrectmaskcoords,resid_toplot);
                axis equal
                shading flat
                title 'Residuals'
                xlabel('xcoords');
                ylabel('ycoords');
                colorbar
            end
        end
        
        
        
        
        
        % Do this through ginput and not roi2poly because the bug prevents
        % roi2poly cursor from clicking at the correct point on the screen.
        function makeTearRegionDisplacementMask(obj)
            obj.makeCustomDisplacementColorPlot();
            disp('Select vertices of the polygon to define the tear region. When complete, click on negative x,y values to terminate the loop.');
            vertex_storage = zeros(0,2);
            while true
                 [x,y] = ginput(1);
                 if x < 0 && y < 0
                     break
                 end
                 vertex_storage(end+1,:) = [x,y];
            end
            
            vertex_storage_pixels = vertex_storage/obj.scan_stepsize + 1;
            vertex_storage_pixels(end+1,:) = vertex_storage_pixels(1,:);
            mask = poly2mask(vertex_storage_pixels(:,1),vertex_storage_pixels(:,2),obj.datacube_size(4),obj.datacube_size(3));
            hold on;
            h = imagesc(obj.xaxis,obj.yaxis,mask);
            h.AlphaData = 0.5;
            obj.tear_mask = mask;
        end
        
        
        
        
        
        % optional SP tolerance should be left empty if not present.
        function detectSolitonWalls(obj,filteramp,filterrange,hotspot_suppression,AA_radius_multiplier,optional_SP_tolerance,optional_AAcircle_min_radius,include_untrimmed)
            if nargin < 2
                filteramp = 0.3;
                filterrange = 9;
            end
            if nargin < 6
                optional_SP_tolerance = [];
            end
            if nargin < 7
                optional_AAcircle_min_radius = [];
            end
            if (nargin < 8), include_untrimmed = false; end;
            if isempty(filteramp) || filterrange == 0 && isempty(hotspot_suppression)
                displacement_field = obj.DSC_fit_storage;
            else
                displacement_field = obj.displacementContinuityFilter(filteramp,filterrange);
            end
            if hotspot_suppression
                displacement_field = obj.suppressHotspots(displacement_field);
            end
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            if nargin < 5
                radius_multiplier = 3;
            else
                radius_multiplier = AA_radius_multiplier;
            end
            AA_lininds_mask = obj.getAACircleMask(radius_multiplier,optional_AAcircle_min_radius,include_untrimmed);
            [fighn] = obj.makeCustomDisplacementColorPlot(SP_hsv,AB_hsv,filteramp,filterrange,1,0,hotspot_suppression);
            plot_flag = 1;
            upsample_flag = 0;
            SP_id = 1;
            [ SP_line1 ] = registerUniqueSPLines( displacement_field, AA_lininds_mask, SP_id, plot_flag, fighn, upsample_flag, optional_SP_tolerance, obj.tear_mask );
            SP_id = 2;
            [ SP_line2 ] = registerUniqueSPLines( displacement_field, AA_lininds_mask, SP_id, plot_flag, fighn, upsample_flag, optional_SP_tolerance, obj.tear_mask );
            SP_id = 3;
            [ SP_line3 ] = registerUniqueSPLines( displacement_field, AA_lininds_mask, SP_id, plot_flag, fighn, upsample_flag, optional_SP_tolerance, obj.tear_mask );
            obj.soliton_walls_unique = {sparse(logical(SP_line1)),sparse(logical(SP_line2)),sparse(logical(SP_line3))};
            obj.soliton_walls_merged = sparse(logical(SP_line1 | SP_line2 | SP_line3));
            % Get the transition vectors for crossing an SP
            % Recall that the flipping ambiguity means we won't know until
            % runtime at any particular AB region whether we should use the
            % positive or negative vector.
            [ v1, v2, ~ ] = getDSCBasisVectors();
            obj.SP_transition_vectors = zeros(3,2);
            obj.SP_transition_vectors(1,:) = v1;
            obj.SP_transition_vectors(2,:) = v1 - v2;
            obj.SP_transition_vectors(3,:) = v2;
        end
        
        
        
        
        
        % To be run after AA and soliton detection
        % AA_inflation_factor is the fraction of the circle fit radius
        function detectABregions(obj,AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin,tear_boundary_trim_margin)
            obj.AB_area_storage = cell(0,1);
            if nargin < 2
                AA_inflation_factor = 0;
            end
            if (nargin < 5), tear_boundary_trim_margin = []; end;
            allmask = true(obj.datacube_size(3:4));
            % Segment out the SPs
            allmask(obj.soliton_walls_merged) = false;
            xbase = 1:obj.datacube_size(4);
            ybase = 1:obj.datacube_size(3);
            [x,y] = meshgrid(xbase,ybase);
            if AA_inflation_factor > 0
                centers = obj.untrimmed_AA_centers;
                if isempty(obj.AA_circlefit_values)
                    radii = obj.AA_Gaussian_Circle_Fit_Params(:,3);
                else
                    radii = obj.AA_circlefit_values(:,3);
                end
                r = mean(radii);
                for i = 1:size(centers,1)
                    % Only inflate false ball if we are close to the edge
                    % of the image where soliton registration tends to
                    % screw up
                    tc = centers(i,:);
                    bm = boundary_inflation_margin;
                    inmargin = tc(1) <= bm | tc(1) > obj.datacube_size(3)-bm | tc(2) <= bm | tc(2) > obj.datacube_size(4)-bm;
                    if inmargin
                        allmask(isInCircle(x,y,tc(1),tc(2),r*AA_inflation_factor)) = false;
                    end
                end
            end
            B = false(size(allmask));
            if boundary_trim_margin > 0
                [ allmask ] = trimArray( allmask,boundary_trim_margin,B );
            end
            if ~isempty(tear_boundary_trim_margin) && ~isempty(obj.tear_mask)
                tear_bmask = boundarymask(obj.tear_mask);
                tear_bmask_wide = tear_bmask;
                for i = 1:tear_boundary_trim_margin
                    tear_bmask_wide = boundarymask(tear_bmask_wide) | tear_bmask_wide;
                    tear_bmask_wide(obj.tear_mask) = false;
                end
                allmask(tear_bmask_wide) = false;
            end
            % NPK removed this portion on 06/26/2020 for better analysis
            % near edges of dataset
            toosmall = 5;
            allmask = bwareaopen(allmask,toosmall,4);
            figure; imagesc(allmask); title('Segmentation mask for AB connectivity detection');
            set(gca,'yDir','normal');
            
            if ~isempty(obj.tear_mask)
                allmask(obj.tear_mask) = false;
            end
            segres = bwconncomp(allmask,4);
            nAA = segres.NumObjects;
            obj.AB_centroids = zeros(nAA,2);
            obj.AB_segmentation_struct = segres;
            obj.AB_label_matrix = bwlabel(allmask,4);
            for i = 1:nAA
                newB = B;
                newB(segres.PixelIdxList{i}) = 1;
                thisAB = sparse(logical(newB));
                obj.AB_area_storage{i} = thisAB;
                obj.AB_centroids(i,:) = mean([x(thisAB), y(thisAB)]);
            end
            obj.makeCustomDisplacementColorPlot([],[],[],[],1);
            hold on
            scatter(obj.AB_centroids(:,1),obj.AB_centroids(:,2),200,'filled','k');
            title('AB centroid positions');
        end
        
        
        
        
        
        function plotGeometryFit(obj,filteramp,filterrange,hotspot_suppression)
            % Primary colors may make this easier to visualize.
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            [figh] = obj.makeCustomDisplacementColorPlot(SP_hsv,AB_hsv,filteramp,filterrange,1,0,hotspot_suppression);
            hold on
            h = imagesc(obj.xaxis,obj.yaxis,full(obj.soliton_walls_merged));
            h.AlphaData = 0.5;
            h2 = imagesc(obj.xaxis,obj.yaxis,label2rgb(obj.AB_label_matrix));
            h2.AlphaData = 0.5;
            scatter(obj.untrimmed_AA_centers(:,1)*obj.scan_stepsize,obj.untrimmed_AA_centers(:,2)*obj.scan_stepsize,100,'filled','r');
            scatter(obj.AB_centroids(:,1)*obj.scan_stepsize,obj.AB_centroids(:,2)*obj.scan_stepsize,100,'filled','k');
            if isempty(obj.AA_circlefit_values)
                AApos = obj.AA_Gaussian_Circle_Fit_Params(:,1:2);
                AArad = obj.AA_Gaussian_Circle_Fit_Params(:,3);
            else
                AApos = obj.AA_circlefit_values(:,1:2);
                AArad = obj.AA_circlefit_values(:,3);
            end
            viscircles(AApos*obj.scan_stepsize,AArad*obj.scan_stepsize,'Color','r');
            
            figure
            imagesc(full(obj.soliton_walls_merged));
            hold on
            h2 = imagesc(label2rgb(obj.AB_label_matrix));
            h2.AlphaData = 0.5;
            scatter(obj.untrimmed_AA_centers(:,1),obj.untrimmed_AA_centers(:,2),100,'filled','r');
            scatter(obj.AB_centroids(:,1),obj.AB_centroids(:,2),100,'filled','k');
            viscircles(AApos,AArad,'Color','r');
            set(gca,'yDir','normal');
        end
        
        
        
        
        
        % Wrapper function.
        function fitAllGeometry(obj,filteramp,filterrange,hotspot_suppression,AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin)
            obj.fitAAToCircle(0,filteramp,filterrange,hotspot_suppression);
            obj.detectSolitonWalls(filteramp,filterrange,hotspot_suppression,AA_inflation_factor);
            obj.detectABregions(AA_inflation_factor,boundary_inflation_margin,boundary_trim_margin);
            obj.plotGeometryFit(filteramp,filterrange,hotspot_suppression);
        end
        
        
        
        
        
        function adjustSolitonWalls(obj)
            obj.plotGeometryFit(100,1,0);
            close
            f = gcf;
%             num = input('Please select a soliton wall (red = 1, green = 2, blue = 3');
            disp('Please click on pixels that should be included in the soliton wall.');
            disp('After each pixel, rezoom if needed and press any key to continue.');
            disp('ginput on negative (x,y) to terminate.');
            disp('These pixels will be added to the merged soliton walls, but not set to any one soliton wall type.');
            while true
                input('Press any key to continue');
                [x,y] = ginput(1);
                if x < 0 && y < 0
                    break 
                end
                cind = round(x/obj.scan_stepsize + 1);
                rind = round(y/obj.scan_stepsize + 1);
                obj.soliton_walls_merged(rind,cind) = true;
            end
            obj.plotGeometryFit(100,1,0);
            disp('The new merged soliton wall registration is displayed.');
        end
        
        
        
        % no default setting on the latter two arguments because you need
        % to be paying attention to not get wrong results.
        function [pixel_sizes,nm_sizes,angle_storage,side_length_storage] = getABareaStatistics(obj,disallow_edge_intersect,boundary_edge_margin)
            nAB = numel(obj.AB_area_storage);
            if ~disallow_edge_intersect
                pixel_sizes = zeros(nAB,1);
                for i = 1:nAB
                    pixel_sizes(i) = nnz(obj.AB_area_storage{i});
                end
                nm_sizes = pixel_sizes*(obj.scan_stepsize^2);  % NPK fixed a very erroneous formula here on 04/28/2020.
                angle_storage = [];
                side_length_storage = [];
            else  % Want to exclude AB domains that intersect with some boundary edge so we can get accurate domain statistics.
                disp('hi');
                falsemat = false(obj.datacube_size(3:4));
                edgemask = trimArray(falsemat,boundary_edge_margin,true(obj.datacube_size(3:4)));
                % cycle through each AB domain and see if it is connected
                % to the boundary of the image. 
                angle_storage = zeros(nAB,3);
                side_length_storage = zeros(nAB,3);
                nm_sizes = zeros(nAB,1);
                pixel_sizes = zeros(nAB,1);
                assigned = false(nAB,1);
                for i = 1:nAB
                    fprintf('Computing statistics on AB domain %d of %d.\n',i,nAB);
                    this_AB = full(obj.AB_area_storage{i});
                    union = this_AB | edgemask;
                    conres = bwconncomp(union);
                    if conres.NumObjects == 2  % then they are disjoint, proceed
                        assigned(i) = true;
                        n_pixels = nnz(this_AB);
                        nm_size = n_pixels*(obj.scan_stepsize^2);
                        nm_sizes(i) = nm_size;
                        pixel_sizes(i) = n_pixels;
                        % Figure out the angles by using information from
                        % corresponding AA vertices. Hard to refactor this
                        % because it won't work for all domains, only those
                        % not bordering edge
                        nAA = size(obj.untrimmed_AA_centers,1);
                        AAdists = [zeros(nAA,1),(1:nAA)'];
                        for j = 1:nAA
                            thisAA = falsemat;
                            % These need to be flipped in order.
                            thisAA(obj.untrimmed_AA_centers(j,2),obj.untrimmed_AA_centers(j,1)) = true;  % single pixel
%                             comptest1 = bwconncomp(this_AB|thisAA);
%                             if comptest1.NumObjects == 1
%                                 AAdists(j) = 0;
%                             else
                                [ ~,distance ] = bwshortestpath( this_AB, thisAA, true, false );
%                                 AAdists(j) = nnz(path_pixels);  % there is always the diagonal issue here, but the correct AA vertices should be so dissimilar from all else that this should be really obvious.
                                AAdists(j) = distance;
%                             end
                        end
                        sAAdists = sortrows(AAdists);
                        % extract the three vertices as those of the
                        % smallest path distances
                        vA_idx = sAAdists(1,2);
                        vB_idx = sAAdists(2,2);
                        vC_idx = sAAdists(3,2);
                        % get pixel [x,y] coordinates and shift to nm.
                        % I think these are in xy, but it doesn't really
                        % matter
                        vA = obj.untrimmed_AA_centers(vA_idx,:)*obj.scan_stepsize;
                        vB = obj.untrimmed_AA_centers(vB_idx,:)*obj.scan_stepsize;
                        vC = obj.untrimmed_AA_centers(vC_idx,:)*obj.scan_stepsize;
                        abdist = sum((vB-vA).^2,2).^0.5;
                        bcdist = sum((vC-vB).^2,2).^0.5;
                        acdist = sum((vC-vA).^2,2).^0.5;
                        Cangle = acos((bcdist^2+acdist^2-abdist^2)/(2*bcdist*acdist));
                        Bangle = acos((bcdist^2+abdist^2-acdist^2)/(2*bcdist*abdist));
                        Aangle = acos((acdist^2+abdist^2-bcdist^2)/(2*acdist*abdist));
                        % since the ordering of the vertices is arbitrary,
                        % store via magnitude
                        these_angles = [Aangle,Bangle,Cangle];
                        these_distances = [abdist,bcdist,acdist];
                        angle_storage(i,:) = rad2deg(sort(these_angles));  % store in degrees
                        side_length_storage(i,:) = sort(these_distances);
                    elseif conres.NumObjects ~= 1
                        error('Unexpected number of connectivity components in getABareaStatistics(). Something has gone wrong.');
                    end  % else do nothing because the AB domain intersects the boundary region
                end
                angle_storage(~assigned,:) = [];
                side_length_storage(~assigned,:) = [];
                nm_sizes(~assigned,:) = [];
                pixel_sizes = [];
                % assign as instance variables in case we want to retrieve
                % later.
                obj.ABangles = angle_storage;
                obj.ABsidelengths = side_length_storage;
                obj.ABareas_nm = nm_sizes;
                obj.ABcalc_boundary_edge_margin = boundary_edge_margin;
                % This last entry serves as a memory of how we got to the
                % stored values.
                % Make plots
                
            end
        end
        
        
        
        
        
        
        % all return values are either in nm or degrees
        function [mean_area,stdev_area,mean_sidelength,stdev_sidelength,...
                mean_smallest_angle,stdev_smallest_angle,mean_angle_difference,...
                stdev_angle_difference] = plotABareaStatistics(obj)
            mean_area = mean(obj.ABareas_nm(:));
            stdev_area = std(obj.ABareas_nm(:));
            mean_sidelength = mean(obj.ABsidelengths(:));
            stdev_sidelength = std(obj.ABsidelengths(:));
            mean_smallest_angle = mean(obj.ABangles(:,1));
            stdev_smallest_angle = std(obj.ABangles(:,1));
            mean_angle_difference = mean(obj.ABangles(:,3)-obj.ABangles(:,1));
            stdev_angle_difference = std(obj.ABangles(:,3)-obj.ABangles(:,1));
            
            figure; histogram(obj.ABareas_nm(:));
            title('AB skeleton domain areas');
            xlabel('nm area');
            ylabel('Counts');
            figure; histogram(obj.ABsidelengths(:));
            title('AB skeleton domain sidelengths');
            xlabel('Nearest-neighbor AA straight line distance');
            ylabel('Counts');
            figure; histogram(obj.ABangles(:,1));
            title('AB skeleton domain smallest angle');
            xlabel('Smallest domain angle');
            ylabel('Counts');
            figure; histogram(obj.ABangles(:,3)-obj.ABangles(:,1));
            title('AB skeleton domain triangular distortion');
            xlabel('Smallest domain angle minus largest domain angle');
            ylabel('Counts');
        end
        
        
        
        
        
        
        function [pixel_radii,nm_radii,nm_average,nm_std] = getAACirclefitStatistics(obj,make_plots,bin_edges)
            pixel_radii = obj.AA_circlefit_values(:,3);
            nm_radii = obj.AA_circlefit_values(:,3)*obj.scan_stepsize;
            if make_plots
                figure
                if nargin < 3
                    histogram(nm_radii);
                else
                    histogram(nm_radii,bin_edges);
                end
                xlabel('Radius of AA center (nm)');
                ylabel('Counts');
                title('Histogram of circle-fit AA domain geometry');
                set(gca,'FontSize',14);
            end
            nm_average = mean(nm_radii,'omitnan');
            nm_std = std(nm_radii,'omitnan');
        end
        
        
        
        
        
        % Omitnans allow post-processing for outliers, as of 05/23/2020
        function [semimajor_nm_av,semimajor_nm_std,semiminor_nm_av,semiminor_nm_std,angle_av,angle_std] = ...
                getAAellipticGaussianFitStatistics(obj,angle_hist_edges)
            semimajor_nm = obj.AAgaussianfit_FWHM_ellipses(:,2)*obj.scan_stepsize;
            semiminor_nm = obj.AAgaussianfit_FWHM_ellipses(:,3)*obj.scan_stepsize;
            angle = rad2deg(obj.AAgaussianfit_FWHM_ellipses(:,1));
            semimajor_nm_av = mean(semimajor_nm,'omitnan');
            semiminor_nm_av = mean(semiminor_nm,'omitnan');
            angle_av = mean(angle,'omitnan');
            semimajor_nm_std = std(semimajor_nm,'omitnan');
            semiminor_nm_std = std(semiminor_nm,'omitnan');
            angle_std = std(angle,'omitnan');
            
            
                figure
                histogram(semimajor_nm);
                xlabel('Semimajor axis of Gaussian FWHM (nm)');
                ylabel('Counts');
                title('Histogram of elliptic Gaussian-fit AA domain geometry');
                set(gca,'FontSize',14);
                
                figure
                histogram(semiminor_nm);
                xlabel('Semiminor axis of Gaussian FWHM (nm)');
                ylabel('Counts');
                title('Histogram of elliptic Gaussian-fit AA domain geometry');
                set(gca,'FontSize',14);
                
                figure
                if nargin < 2
                    histogram(angle);
                else
                    histogram(angle,angle_hist_edges);
                end
                xlabel('Rotation angle of Gaussian (degrees)');
                ylabel('Counts');
                title('Histogram of elliptic Gaussian-fit AA domain geometry');
                set(gca,'FontSize',14);
            
        end
        
        
        
        
        
        function plotAACircleFitOnCustomColor(obj,filteramp,filterrange,suppress_hotspots,ignore_stepsize)
            if (nargin < 5), ignore_stepsize = false; end;
            obj.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,suppress_hotspots);
            hold on
            nAA = size(obj.AA_circlefit_values,1);
            for i = 1:nAA
                cent = obj.AA_circlefit_values(i,1:2);
                radius = obj.AA_circlefit_values(i,3);
                if ~ignore_stepsize
                    % NPK modified 05/29/2020 to properly adjust for the
                    % axes of the plot, starting at zero.
                    cent = (cent-1)*obj.scan_stepsize;
                    radius = radius*obj.scan_stepsize;
                end
                viscircles(cent,radius,'Color','g');
            end
            axis equal
        end
        
        
        
        
        
        function plotAAGaussianEllipseFitOnCustomColor(obj,filteramp,filterrange,suppress_hotspots,ignore_stepsize)
            if (nargin < 5), ignore_stepsize = false; end;
            obj.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,suppress_hotspots);
            hold on
            nAA = size(obj.AAgaussianfit_FWHM_ellipses,1);
            for i = 1:nAA
                x0 = obj.AAgaussianfit_FWHM_ellipses(i,4);
                y0 = obj.AAgaussianfit_FWHM_ellipses(i,5);
                semimajor = obj.AAgaussianfit_FWHM_ellipses(i,2);
                semiminor = obj.AAgaussianfit_FWHM_ellipses(i,3);
                if ~ignore_stepsize
                    x0 = (x0-1)*obj.scan_stepsize;
                    y0 = (y0-1)*obj.scan_stepsize;
                    semimajor = semimajor*obj.scan_stepsize;
                    semiminor = semiminor*obj.scan_stepsize;
                end
                angle = obj.AAgaussianfit_FWHM_ellipses(i,1);
%                 ellipse(cent,radius,'Color','g');
                h=ellipse(semimajor,semiminor,angle,x0,y0);
                h.Color = 'g';
                h.LineWidth = 1.5;
            end
            axis equal
        end
        
        
        
        % Circumstantially, it appears that each point only works
        % topologically with either the v1 or the v2 basis. So try v1
        % first, show user the results, discard if discontinuous, and then
        % use the v2 starting point.
        function setEmitterOrientations(obj,direction,tear_protocol)
            % While AB centroids are in xy coordinates, we will manipulate
            % the emitters on a pixel level, so fliplr to get into that
            % orientation.
            if (nargin < 3), tear_protocol = false; end;
            obj.emitter_pixel_pos = fliplr(round(obj.AB_centroids));
            nAB = size(obj.AB_centroids,1);
            if ~tear_protocol || isempty(obj.emitter_displacements)
                obj.emitter_displacements = zeros(size(obj.emitter_pixel_pos));
                obj.emitter_is_set = false(nAB,1);
            else
                init_emitter_displacements = obj.emitter_displacements;
                init_isset = obj.emitter_is_set;
            end
            
            
            
            % Original approach: use largest AB region
            [pixel_sizes] = obj.getABareaStatistics(0,[]);
            [~,first_idx] = max(pixel_sizes);
            randcount = 0;
            while true
                if randcount
                    first_idx = randi(size(obj.AB_centroids,1));
                end
                if tear_protocol
                    obj.plotGeometryFit(100,1,0);
                    close
                    disp('Please click on the starting AB emitter. Will snap to centroid.');
                    [xclick,yclick] = ginput(1);
                    rindex = yclick/obj.scan_stepsize+1;
                    cindex = xclick/obj.scan_stepsize+1;
                    rdist = obj.emitter_pixel_pos(:,1) - rindex;
                    cdist = obj.emitter_pixel_pos(:,2) - cindex;
                    dists = sqrt(rdist.^2 + cdist.^2);
                    [~,first_idx] = min(dists);
                end
                randcount = 1;
                % Define an arbitrary AB displacement for the largest AB region.
                [ v1, v2, ~ ] = getDSCBasisVectors();
                %             w1 = v1 + v2;
                %             w2 = 2*v2 - v1;
                if nargin > 1
                    switch direction
                        case 'v1'
                            first_disp = v1;
                        case 'v2'
                            first_disp = v2;
                    end
                else
                    first_disp = v2;
                end
                obj.emitter_displacements(first_idx,:) = first_disp;
                obj.emitter_is_set(first_idx) = true;
                
                % Orders the population of emitter displacements
                obj.updateNeighboringDisplacementsOrderControl(first_idx);
                % Handle class update on all relevant properties
                plotscale = 4;
                figure; quiver(obj.emitter_pixel_pos(:,2),obj.emitter_pixel_pos(:,1),plotscale*obj.emitter_displacements(:,1),plotscale*obj.emitter_displacements(:,2),0);
                title('Extended zone displacement emitters, v1 starting point');
                disp('If this vector field exhibits a discontinuity, then the alternative starting displacement should be chosen.');
                yn = input('Do you wish to erase this calculation and choose the v1 starting point? 1/0 ');
                if yn
                    if tear_protocol
                        obj.emitter_displacements = init_emitter_displacements;
                        obj.emitter_is_set = init_isset;
                    else
                        obj.emitter_displacements = zeros(size(obj.emitter_pixel_pos));
                        obj.emitter_is_set = false(nAB,1);
                    end
                    first_disp = v1;
                    obj.emitter_displacements(first_idx,:) = first_disp;
                    obj.emitter_is_set(first_idx) = true;
                    obj.updateNeighboringDisplacementsOrderControl(first_idx);
                    
%                     figure; quiver(obj.emitter_pixel_pos(:,2),obj.emitter_pixel_pos(:,1),plotscale*obj.emitter_displacements(:,1),plotscale*obj.emitter_displacements(:,2),0);
                    obj.makeCustomDisplacementColorPlot(); hold on;
                    h = imagesc(obj.xaxis,obj.yaxis,full(obj.soliton_walls_merged));
                    h.AlphaData = 0.3;
                    plotscale = 2;
                    quiver((obj.emitter_pixel_pos(:,2)-1)*obj.scan_stepsize,(obj.emitter_pixel_pos(:,1)-1)*obj.scan_stepsize,plotscale*obj.emitter_displacements(:,1),plotscale*obj.emitter_displacements(:,2),0,'Color','w','LineWidth',2);
                    title('Extended zone displacement emitters, v1 starting point');
                    disp('The v1 starting point extended zone field will be saved.');
                end
                yn2 = input('Do you wish to try a different starting index? 1/0');
                if ~yn2
                    break
                end
            end
            disp('Exiting method setEmitterOrientations()');
        end
        
        
        
        
        
        % Prerequisite is running the setEmitterOrientations() method.
        % Can set makePlots to 0 (nothing at all), 1 (quiver plots only),
        % or 2 (quiver and energy plots)
        function [RMSgradient_energy_storage, magnetic_energy_storage, num_vectors_changed_storage] = ...
                annealDisplacementField(obj, filteramp, filterrange, useParallel, makePlots, scheduleID)
            tic
            if isempty(filteramp) || filterrange == 0
                dfield = obj.DSC_fit_storage;
            else
                dfield = obj.displacementContinuityFilter(filteramp,filterrange);
            end
            obj.annealed_dfield = zeros(obj.datacube_size(3),obj.datacube_size(4),2);
            
            % Visualize initial displacement field
            [r,c,h] = size(dfield);
            rbase = 1:r;
            cbase = 1:c;
            [rgrid,cgrid] = meshgrid(rbase,cbase);
            magnets = cell(r,c);
            energies = zeros(r,c);
            fixed_emitters_indices = obj.emitter_pixel_pos;  % This needs to be swtiched because we will use it to compare distances against indices
            fixed_emitters_directions = obj.emitter_displacements;
            figure
            dfieldx = dfield(:,:,1);
            dfieldy = dfield(:,:,2);
            quiver(rgrid(:),cgrid(:),dfieldx(:),dfieldy(:),0);
            title(sprintf('Starting quiver'));
            
            % Set the annealing schedule
            switch scheduleID
                case 1 % Original schedule (mostly) implemented on NPK's laptop, though last iteration failed.
                    niter = 6;
                    f_coupling_schedule = [100000,2;
                        10000,2;
                        100,1.1;
                        10,1.1;
                        1,1.1;
                        0,1.1];
                    n_coupling_schedule = [0;
                        1;
                        10;
                        100;
                        1000;
                        100000];
                case 2
                    niter = 12;
                    f_coupling_schedule = horzcat([logspace(5,0,11)';0],ones(12,1));
                    n_coupling_schedule = [0;logspace(0,5,11)'];
                case 3
                    niter = 12;
                    f_coupling_schedule = horzcat([logspace(5,0,11)';0],3*ones(12,1));
                    n_coupling_schedule = [0;logspace(0,5,11)'];
                case 4
                    niter = 5;
                    f_coupling_schedule = vertcat([100000,2],zeros(4,2));
                    n_coupling_schedule = [0;1;1;1;1];
            end
            
            % Build the cell array of DisplacementMagnets;
            permutation_number = 5;
            PAD = 2;
            xbounds = zeros(1,2);
            ybounds = zeros(1,2);
            xbounds(1) = max(obj.emitter_displacements(:,1)) + PAD;
            xbounds(2) = min(obj.emitter_displacements(:,1)) - PAD;
            ybounds(1) = max(obj.emitter_displacements(:,2)) + PAD;
            ybounds(2) = min(obj.emitter_displacements(:,2)) - PAD;
            disp('Building displacement magnets...');
            ub2 = obj.datacube_size(4);
            parfor i = 1:obj.datacube_size(3)
                for j = 1:ub2
                    magnets{i,j} = DisplacementEquivalenceClassMagnet([dfieldx(i,j),dfieldy(i,j)],[i,j],permutation_number,xbounds,ybounds);
                end
            end
            
            % Compute starting energies
            f_coupling = f_coupling_schedule(1,:);
            n_coupling = n_coupling_schedule(1);
            disp('Computing starting energies...');
            for i = 1:obj.datacube_size(3)
                for j = 1:ub2
                    energies(i,j) = magnets{i,j}.evaluateCurrentEnergy(dfield,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
                end
            end
            figure
            imagesc(energies); colormap(fire);
            starting_energy = sum(sum(energies));
            colorbar;
            title(sprintf('Beginning Energy, energy %f',starting_energy));
            
            RMSgradient_energy_storage = zeros(niter+1,1);
            magnetic_energy_storage = zeros(niter+1,1);
            num_vectors_changed_storage = zeros(niter+1,1);
            [RMSgradient_energy,magnetic_energy,num_vectors_changed] = getAnnealedDfieldStandardEnergy(obj,obj.DSC_fit_storage,magnets);
            RMSgradient_energy_storage(1) = RMSgradient_energy;
            magnetic_energy_storage(1) = magnetic_energy;
            num_vectors_changed_storage(1) = num_vectors_changed;
            
            transition_flag = 1;
            for q = 1:niter
                f_coupling = f_coupling_schedule(q,:);
                n_coupling = n_coupling_schedule(q);
                %                 if q > 1
                %                     f_coupling = [10000,2];
                %                     n_coupling = 1;
                %                 end
                %                 if q > 2
                %                     f_coupling = [100,1.1];
                %                     n_coupling = 10;
                %                 end
                %                 if q > 3
                %                     f_coupling = [10,1.1];
                %                     n_coupling = 100;
                %                 end
                if useParallel
                    fprintf('Transitioning vectors for iteration %d of %d.\n',q,niter);
                    parfor i = 1:obj.datacube_size(3)
                        for j = 1:ub2
                            %                 energies(i,j) = magnets{i,j}.evaluateCurrentEnergy(n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
                            [energies(i,j),~] = magnets{i,j}.transitionToBestNewVector(dfield,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling,transition_flag);
                        end
                    end
                else
                    for i = 1:obj.datacube_size(3)
                        fprintf('Computing flip row %d of %d for iteration %d of %d.\n',i,obj.datacube_size(3),q,niter);
                        for j = 1:obj.datacube_size(4)
                            %                 energies(i,j) = magnets{i,j}.evaluateCurrentEnergy(n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling);
                            [energies(i,j),~] = magnets{i,j}.transitionToBestNewVector(dfield,n_coupling,fixed_emitters_directions,fixed_emitters_indices,f_coupling,transition_flag);
                        end
                    end
                end
                
                if makePlots > 1
                    figure
                    imagesc(energies); colormap(fire);
                    this_energy = sum(sum(energies));
                    hold on;
                    colorbar;
                    title(sprintf('Iteration %d, energy %f',q,this_energy));
                end
                if makePlots > 0
                    figure
                    dfieldx = obj.annealed_dfield(:,:,1);
                    dfieldy = obj.annealed_dfield(:,:,2);
                    quiver(rgrid(:),cgrid(:),dfieldx(:),dfieldy(:),0);
                    hold on
                    quiver(fixed_emitters_indices(:,2),fixed_emitters_indices(:,1),fixed_emitters_directions(:,1),fixed_emitters_directions(:,2),0,'r');
                    title(sprintf('Iteration %d',q));
                end
                
                if q == 1
                    last_dfield = obj.DSC_fit_storage;
                end
                for i = 1:obj.datacube_size(3)
                    for j = 1:obj.datacube_size(4)
                        obj.annealed_dfield(i,j,:) = permute(magnets{i,j}.current_displacement,[1,3,2]);
                    end
                end
                dfield = obj.annealed_dfield;
                
                % Compute standardized energy, number of significant
                % changes from the last displacement field iteration.
                [RMSgradient_energy,magnetic_energy,num_vectors_changed] = obj.getAnnealedDfieldStandardEnergy(last_dfield,magnets);
                RMSgradient_energy_storage(q+1) = RMSgradient_energy;
                magnetic_energy_storage(q+1) = magnetic_energy;
                num_vectors_changed_storage(q+1) = num_vectors_changed;
                last_dfield = dfield;
            end
            finaltime = toc;
            fprintf('Annealing schedule id#%d completed in %.2f minutes with a final gradient RMS energy of %.2f.\n',scheduleID,finaltime,RMSgradient_energy);
        end
        
        
        
        
        
        % 04/26/2020 development function to try to figure out what is not
        % quite working about the original annealing script.
        %
        function annealDisplacementField2(obj, filteramp, filterrange, makePlots)
            
            dfield = obj.displacementContinuityFilter(filteramp,filterrange);
            dfieldx = dfield(:,:,1);
            dfieldy = dfield(:,:,2);
            dfieldxlin = dfieldx(:);
            dfieldylin = dfieldy(:);
            
            perm_val = 10;
            basis = permn(-perm_val:perm_val,2);
            [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
            w1 = v1 + v2;
            w2 = 2*v2 - v1;
            W = [w1',w2'];
            genpoints = W*basis';
            % Trim down the permutation basis on the grounds of emitter
            % locations
            xmax = max(obj.emitter_displacements(:,1));
            xmin = min(obj.emitter_displacements(:,1));
            ymax = max(obj.emitter_displacements(:,2));
            ymin = min(obj.emitter_displacements(:,2));
            padval = 4;
            remove = genpoints(1,:) > xmax+padval | genpoints(1,:) < xmin-padval | genpoints(2,:) > ymax+padval | genpoints(2,:) < ymin-padval;
            genpoints(:,remove) = [];
            genpointsx = genpoints(1,:);
            genpointsy = genpoints(2,:);
            
            genpointsx_plus = genpointsx + dfieldxlin;
            genpointsx_minus = genpointsx - dfieldxlin;
            equivalent_x_disps = [genpointsx_plus,genpointsx_minus,-dfieldxlin];
            
            genpointsy_plus = genpointsy + dfieldylin;
            genpointsy_minus = genpointsy - dfieldylin;
            equivalent_y_disps = [genpointsy_plus,genpointsy_minus,-dfieldylin];
            
            % Make this plot to ensure that the permutation basis has
            % sufficiently covered the emitter range.
            figure; scatter(equivalent_x_disps(2,:),equivalent_y_disps(2,:)); hold on; plotFullDisplacementHexagons(gca);
            hold on
            scatter(obj.emitter_displacements(:,1),obj.emitter_displacements(:,2),'r');
            scatter(genpointsx,genpointsy,'g');
            
            % set initial displacements to the emitter of that segmented AB
            % region.
            initial_xdisp = zeros(size(equivalent_x_disps,1),1);
            initial_ydisp = zeros(size(equivalent_y_disps,1),1);
            nAB = obj.AB_segmentation_struct.NumObjects;
            fprintf('Guiding initial displacements to emitters...\n');
            for i = 1:nAB
                this_emitterdisp = obj.emitter_displacements(i,:);
                lininds = obj.AB_segmentation_struct.PixelIdxList{i};
                % This works because row dimension is the linear
                % representation of the 2D mat into which segmentation
                % assignments were made.
                useable_xdisps = equivalent_x_disps(lininds,:);
                useable_ydisps = equivalent_y_disps(lininds,:);
                dists = ((useable_xdisps-this_emitterdisp(1)).^2 + (useable_ydisps-this_emitterdisp(2)).^2).^0.5;
                [~,idxs] = min(dists,[],2);
                [r,c] = size(useable_ydisps);
                linidxs = sub2ind([r,c],(1:r)',idxs);
                initial_xdisp(lininds,:) = useable_xdisps(linidxs);
                initial_ydisp(lininds,:) = useable_ydisps(linidxs);
            end
            
            % If not assigned to an emitter, compute distance to emitters,
            % and then go with the closest emitter
            not_assigned = initial_xdisp == 0 & initial_ydisp == 0;
            useable_xdisps = equivalent_x_disps(not_assigned,:);
            useable_ydisps = equivalent_y_disps(not_assigned,:);
            [rf,cf,hf] = size(obj.DSC_fit_storage);
            % NPK 05/24/2020 changed from [cb,rb] = meshgrid(1:rf,1:cf);
            % There seemed to be a dimension problem with dataset 10
%             [rb,cb] = meshgrid(1:rf,1:cf);  % Try 2. This didn't work
            [cb,rb] = meshgrid(1:cf,1:rf);  % Try 3. This seems correct. Works on DS 10. Also works on DS 18, with equivalent results.
            rblin = rb(:);
            cblin = cb(:);
            %             pixeldists = ((rblin(not_assigned)-obj.emitter_pixel_pos(:,1)').^2+(cblin(not_assigned)-obj.emitter_pixel_pos(:,2)').^2).^0.5;
            %             pixeldists = ((rblin(not_assigned)-rblin(~not_assigned)').^2+(cblin(not_assigned)-cblin(not_assigned)').^2).^0.5;
            %             [~,emitteridx] = min(pixeldists,[],2);
            %             rnlin = rblin(~not_assigned);
            %             cnlin = cblin(~not_assigned);
            %             rnlin(
            %             xemitterdisp = obj.emitter_displacements(emitteridx,1);
            %             yemitterdisp = obj.emitter_displacements(emitteridx,2);
            % Now that we have the desired displacements, evaluate those
            % against the equivalent displacements for each point.
            % Vectorized over unassigned pixels
            % %
            % %             dispdists = ((useable_xdisps-xemitterdisp).^2+(useable_ydisps-yemitterdisp).^2).^0.5;
            % %             [~,idxs] = min(dispdists,[],2);
            % %             [r,c] = size(useable_ydisps);
            % %             linidxs = sub2ind([r,c],(1:r)',idxs);
            % %             initial_xdisp(not_assigned) = useable_xdisps(linidxs);
            % %             initial_ydisp(not_assigned) = useable_ydisps(linidxs);
            % Get neighboring linear indices
            [nbaser,nbasec] = meshgrid(-1:1,-1:1);
            cn = cblin + nbasec(:)';
            rn = rblin + nbaser(:)';
            delete = cn < 1 | cn > cf | rn < 1 | rn > rf;
            cn(delete) = nan;
            rn(delete) = nan;
            neighbor_lininds = sub2ind([rf,cf],rn,cn);
            % This loop will replace nans with neighbors.
            for k = 1:8
                [rnanres,cnanres] = find(isnan(neighbor_lininds));
                cnanres_shifted = mod(cnanres,size(neighbor_lininds,2))+1;
                shifted_inds = sub2ind(size(neighbor_lininds),rnanres,cnanres_shifted);
                neighbor_lininds(isnan(neighbor_lininds)) = neighbor_lininds(shifted_inds);
            end
            
            
            is_neighbor_assigned_mat = ~not_assigned(neighbor_lininds);
            xneighbors = initial_xdisp(neighbor_lininds);
            yneighbors = initial_ydisp(neighbor_lininds);
            xneighbors(~is_neighbor_assigned_mat) = nan;
            yneighbors(~is_neighbor_assigned_mat) = nan;
            xmeanndisp = mean(xneighbors,2,'omitnan');
            ymeanndisp = mean(yneighbors,2,'omitnan');
            
            not_assigned_now = not_assigned;
            while nnz(not_assigned_now) > 0
                not_assigned_now = isnan(xmeanndisp);
                is_neighbor_assigned_now_mat = ~not_assigned_now(neighbor_lininds);
                xneighbors = xmeanndisp(neighbor_lininds);
                yneighbors = ymeanndisp(neighbor_lininds);
                xneighbors(~is_neighbor_assigned_now_mat) = nan;
                yneighbors(~is_neighbor_assigned_now_mat) = nan;
                xmeanndisp = mean(xneighbors,2,'omitnan');
                ymeanndisp = mean(yneighbors,2,'omitnan');
            end
            %             ymeanndisp = mean(yneighbors,2,'omitnan');
            % Make any possible assignments
            %             xmeanndispmat = reshape(xmeanndisp,rf,cf);
            %             figure; imagesc(xmeanndispmat); title 'x mean neighbor disp'; colorbar; set(gca,'yDir','normal'); axis equal;
            
            % For all those lacking deterministic assignment, make it on
            % the basis of these extended means as emitters.
            xemittertarget = xmeanndisp(not_assigned);
            yemittertarget = ymeanndisp(not_assigned);
            dists2 = ((useable_xdisps - xemittertarget).^2 + (useable_ydisps - yemittertarget).^2).^0.5;
            [~,emitteridx2] = min(dists2,[],2);
            [r,c] = size(useable_xdisps);
            linidxs = sub2ind([r,c],(1:r)',emitteridx2);
            
            %             dispdists2 = ((useable_xdisps-xemitterdisp).^2+(useable_ydisps-yemitterdisp).^2).^0.5;
            % %             [~,idxs] = min(dispdists,[],2);
            % %             [r,c] = size(useable_ydisps);
            % %             linidxs = sub2ind([r,c],(1:r)',idxs);
            % %             initial_xdisp(not_assigned) = useable_xdisps(linidxs);
            % %             initial_ydisp(not_assigned) = useable_ydisps(linidxs);
            
            initial_xdisp(not_assigned) = useable_xdisps(linidxs);
            initial_ydisp(not_assigned) = useable_ydisps(linidxs);
            if makePlots > 0
                unwrapAndPlotDfield('starting assignments',initial_xdisp,initial_ydisp,[rf,cf]);  % org rf,cf
            end
            %             figure;
            %             scatter(initial_xdisp,initial_ydisp,5,'filled');
            %             initial_xdisps_mat = reshape(initial_xdisp,rf,cf);
            %             initial_ydisps_mat = reshape(initial_ydisp,rf,cf);
            %             figure;
            %             imagesc(initial_xdisps_mat); title 'x displacement'; colorbar; set(gca,'yDir','normal'); axis equal;
            %             figure; imagesc(initial_ydisps_mat); title 'y displacement'; colorbar; set(gca,'yDir','normal'); axis equal;
            
            
            % This is looking very good. Now, compute the ten nearest
            % displacement vectors. We will need these for any sort of a
            % flipping algorithm.
            %
            % Remember too that we already have neighbor lininds! Can use
            % those for the similarity calculation.
            %
            % At the end of this, verify that we have not accidentally
            % flipped to any inequivalent vectors.
            dist_to_chosen = ((equivalent_x_disps - initial_xdisp).^2 + (equivalent_y_disps - initial_ydisp).^2).^0.5;
            [~,I] = sort(dist_to_chosen,2,'ascend');
            rinds = (1:size(I,1))'.*ones(size(I,1),10);
            cinds = I(:,1:10);
            linds = sub2ind(size(I),rinds,cinds);
            reduced_equivalent_x_disps = equivalent_x_disps(linds);
            reduced_equivalent_y_disps = equivalent_y_disps(linds);
            % Okay, these should have the ten closest vectors for any pixel
            % coordinate. So we will use these for any sort of flipping
            % algorithms. Looks correct
            
            % Try deterministic rounds and see how far this gets us:
            %             npixels = size(initial_xdisp,1);
            niter = 10;
            neighbor_lininds_perm = permute(neighbor_lininds,[1,3,2]);
            current_xdisps = initial_xdisp;
            current_ydisps = initial_ydisp;
            num_changed_stor = zeros(niter,1);
            for i = 1:niter
                fprintf('Annealing iteration %d of %d.\n',i,niter);
                current_xdisps_perm = current_xdisps(neighbor_lininds_perm);
                current_ydisps_perm = current_ydisps(neighbor_lininds_perm);
                dists = ((current_xdisps_perm - reduced_equivalent_x_disps).^2 + (current_ydisps_perm - reduced_equivalent_y_disps).^2).^0.5;
                sumdist = sum(dists,3);
                [~,bestcol] = min(sumdist,[],2);
                [rred,cred] = size(reduced_equivalent_x_disps);
                linds = sub2ind([rred,cred],(1:rred)',bestcol);
                changed = current_xdisps ~= reduced_equivalent_x_disps(linds) | current_ydisps ~= reduced_equivalent_y_disps(linds);
                num_changed_stor(i) = nnz(changed);
                current_xdisps = reduced_equivalent_x_disps(linds);
                current_ydisps = reduced_equivalent_y_disps(linds);
                if makePlots > 1 || (makePlots > 0 && i == niter)
                    unwrapAndPlotDfield(sprintf('iteration %d',i),current_xdisps,current_ydisps,[rf,cf]);
                end
            end
            num_changed_stor
            current_xdisps_wrapped = reshape(current_xdisps,rf,cf);
            current_ydisps_wrapped = reshape(current_ydisps,rf,cf);
            if ~isempty(obj.tear_mask)
                current_xdisps_wrapped(obj.tear_mask) = nan;
                current_ydisps_wrapped(obj.tear_mask) = nan;
            end
            obj.annealed_dfield = cat(3,current_xdisps_wrapped,current_ydisps_wrapped);
            
            
            function unwrapAndPlotDfield(titlestring,xdfieldlin,ydfieldlin,siz)
                figure;
                scatter(xdfieldlin,ydfieldlin,5,'filled');
                current_xdisps_mat = reshape(xdfieldlin,siz(1),siz(2));
                current_ydisps_mat = reshape(ydfieldlin,siz(1),siz(2));
                figure; imagesc(current_xdisps_mat); title(sprintf('X displacement %s',titlestring)); colorbar; set(gca,'yDir','normal'); axis equal;
                figure; imagesc(current_ydisps_mat); title(sprintf('Y displacement %s',titlestring)); colorbar; set(gca,'yDir','normal'); axis equal;
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        % Helper function for annealing
        function [RMSgradient_energy,magnetic_energy,num_vectors_changed] = getAnnealedDfieldStandardEnergy(obj,last_dfield,magnetic_cell)
            TOL = 1e-4;
            [gxx,gxy] = gradient(obj.annealed_dfield(:,:,1));
            [gyx,gyy] = gradient(obj.annealed_dfield(:,:,2));
            RMSgradient_energy = sqrt(sum(sum(gxx.^2 + gxy.^2 + gyx.^2 + gyy.^2))/(4*numel(gxx)));
            if nargin > 2
                fprintf('Computing standardized magnetic energies for n_coupling = 1...\n');
                energies = zeros(obj.datacube_size(3),obj.datacube_size(4));
                for i = 1:obj.datacube_size(3)
                    for j = 1:obj.datacube_size(4)
                        energies(i,j) = magnetic_cell{i,j}.evaluateCurrentEnergy(obj.annealed_dfield,1,[],[],[0,0]);
                    end
                end
                magnetic_energy = rms(rms(energies));
            else
                magnetic_energy = [];
            end
            if nargin > 1
                xchange = abs(obj.annealed_dfield(:,:,1) - last_dfield(:,:,1)) > TOL;
                ychange = abs(obj.annealed_dfield(:,:,2) - last_dfield(:,:,2)) > TOL;
                totalchange = xchange | ychange;
                num_vectors_changed = nnz(totalchange);
            else
                num_vectors_changed = [];
            end
        end
        
        
        
        
        function plotAnnealedDisplacementField(obj)
            figure
            imagesc(obj.annealed_dfield(:,:,1));
            title('Annealed x displacement field');
            colormap(parula)
            colorbar;
            axis equal
            set(gca,'yDir','normal');
            figure
            imagesc(obj.annealed_dfield(:,:,2));
            title('Annealed y displacement field');
            colormap(parula)
            colorbar;
            axis equal
            set(gca,'yDir','normal');
        end
        
        
        
        
        
        
        function cropped_annealed_dfield = getCroppedAnnealedDfield(obj,crop_margin,make_plots,rrange,crange)
            if (nargin < 3), make_plots = true; end;
            if (nargin < 4), rrange = []; crange = []; end;
            if isempty(rrange)
                cropped_annealed_dfield(:,:,1) = trimArray(obj.annealed_dfield(:,:,1),crop_margin);
                cropped_annealed_dfield(:,:,2) = trimArray(obj.annealed_dfield(:,:,2),crop_margin);
            else  % xrange, yrange should be a two-element vector each, holding min and max indices in 1 and 2
                [r,c,h] = size(obj.annealed_dfield);
                if rrange(1) < 1
                    rrange(1) = 1;
                end
                if rrange(2) > r
                    rrange(2) = r;
                end
                if crange(1) < 1
                    crange(1) = 1;
                end
                if crange(2) > c
                    crange(2) = c;
                end
                cropped_annealed_dfield = obj.annealed_dfield(rrange(1):rrange(2),crange(1):crange(2),:);
            end
            if make_plots
                figure;
                imagesc(cropped_annealed_dfield(:,:,1));
                title('x dfield cropped');
                set(gca,'yDir','normal');
                colorbar;
                figure;
                imagesc(cropped_annealed_dfield(:,:,2));
                title('y dfield cropped');
                set(gca,'yDir','normal');
                colorbar;
            end
        end
        
        
        
        
        
        
        function twist_angle = fitMoireAngle(obj,crop_margin,rrange,crange,make_plots)
            if nargin < 3
                rrange = [];
                crange = [];
            end
            if nargin < 5
                make_plots = true;
            end
            cropped_annealed_dfield = obj.getCroppedAnnealedDfield(crop_margin,make_plots,rrange,crange);
            cadx = cropped_annealed_dfield(:,:,1);
            cady = cropped_annealed_dfield(:,:,2);
            uX = cadx(:)/10;  % need to convert from angstrom to nm
            uY = cady(:)/10;
            uAMPsq = uX.^2 + uY.^2;
            obj.reformAxes();  % as a precaution, in case they are incorrectly set
            [x0s_in,y0s_in] = meshgrid(obj.xaxis,obj.yaxis);
            if isempty(rrange)
                x0s_in = trimArray(x0s_in,crop_margin);
                y0s_in = trimArray(y0s_in,crop_margin);
            else
                x0s_in = x0s_in(rrange(1):rrange(2),crange(1):crange(2));
                y0s_in = y0s_in(rrange(1):rrange(2),crange(1):crange(2));
            end
            [rsize,csize] = size(x0s_in);
            x0s_in = x0s_in(:);
            y0s_in = y0s_in(:);
%             fitfun = @(thetam) moireFitFun(uX,uY,x0s_in,y0s_in,thetam);
%             fitfun2 = @(params) moireFitFun2(uX,uY,x0s_in,y0s_in,params(1),params(2),params(3));
%             myfun3 = @(xor,yor) rms(moireFitFun2(uX,uY,x0s_in,y0s_in,xor,yor,1.2));
            fitfun3 = @(params) moireFitFun3(uAMPsq,x0s_in,y0s_in,params(1),params(2),params(3));
%             xbase = -200:10:200;
%             ybase = -200:10:200;
%             [xspace,yspace] = meshgrid(xbase,ybase);
%             for i = 1:size(xspace,1)
%                 for j = 1:size(xspace,2)
%                     rmsrsurface(i,j) = myfun3(xspace(i,j),yspace(i,j));
%                 end
%             end
%             figure;
%             surf(xbase,ybase,rmsrsurface);
            
            % multistart needed?
            %%%%%%%%%%%%%%%%%%%%%%% NPK cahnge 05/14/2020 to remove
            %%%%%%%%%%%%%%%%%%%%%%% multistart, consider putting back
            %%%%%%%%%%%%%%%%%%%%%%% later.
%             init_guesses = 0.2:0.2:1.6;
            init_guesses = 1.2;
            thetam_calc_stor = zeros(numel(init_guesses),3);
%             thetam_calc_stor = zeros(numel(init_guesses),1);
            resid_stor = zeros(numel(init_guesses),1);
            origin_guess = [0,120];
            options = optimset('lsqnonlin');
            options.MaxFunEvals = 5000;
            options.Display = 'off';
            for i = 1:numel(init_guesses)
                thetam0 = init_guesses(i);
                %                 init_guess = [origin_guess,thetam0];
                %                 params_calc = lsqnonlin(fitfun2,init_guess);
                init_guess = [origin_guess,thetam0];
                params_calc = lsqnonlin(fitfun3,init_guess,[],[],options);
                %                 thetam_calc = lsqnonlin(fitfun,thetam0);
                %                 thetam_calc_stor(i,:) = params_calc;
                %                 residval = rms(fitfun2(params_calc));
%                 thetam_calc_stor(i) = thetam_calc;
                thetam_calc_stor(i,:) = params_calc;
%                 residval = rms(fitfun(thetam_calc));
                residval = rms(fitfun3(params_calc));
                resid_stor(i) = residval;
            end
%             thetam_calc_stor
%             resid_stor
            
            twist_angle = median(thetam_calc_stor(:,3));
            
            pred_disp_field = moirePredFun3(x0s_in,y0s_in,params_calc(1,1),params_calc(1,2),params_calc(1,3));
            cropped_amp = sum(cropped_annealed_dfield.^2,3).^0.5;
            if isempty(rrange)
                pred_disp_field = reshape(pred_disp_field,[rsize,csize]);
            else
                pred_disp_field = reshape(pred_disp_field,[numel(rrange(1):rrange(2)),numel(crange(1):crange(2))]);
            end
            new_adisp_field = cropped_amp - 10*pred_disp_field.^0.5;
            
%             mean(mean(new_adisp_field))
            
            if make_plots
                figure;
                imagesc(new_adisp_field);
                title('Residual displacement field after Moire fit');
                set(gca,'yDir','normal');
                axis equal
                cbh = colorbar;
                colormap(evalin('base','cmap'));
                ylabel(cbh,'Nonlinear displacement amplitude');
            end
        end
        
        
        
        
        
        function [twist_storage,rrange_storage,crange_storage] = fitMovingMoireAngle(obj,rrange_global,crange_global,window_size,moving_inc,make_plots)
            if nargin < 6
                make_plots = false;
            end
            % try 40 x 40 in 10 pixel movements.
%             [20,160],[40,180]
            
%             rranges = [(20:10:120)',(60:10:160)'];
%             cranges = [(40:10:130)',(80:10:170)'];
            rranges_lower = rrange_global(1):moving_inc:(rrange_global(2)-window_size+1);
            cranges_lower = crange_global(1):moving_inc:(crange_global(2)-window_size+1);
            rranges_upper = rranges_lower+window_size-1;
            cranges_upper = cranges_lower+window_size-1;
            rranges = [rranges_lower(:),rranges_upper(:)];
            cranges = [cranges_lower(:),cranges_upper(:)];
            nr = size(rranges,1);
            nc = size(cranges,1);
            twist_storage = zeros(nr,nc);
            rrange_storage = zeros(nr,nc,2);
            crange_storage = zeros(nr,nc,2);
            for i = 1:nr
                fprintf('Fitting row %d of %d.\n',i,nr);
                for j = 1:nc
                    this_rrange = rranges(i,:);
                    this_crange = cranges(j,:);
                    this_twist = obj.fitMoireAngle([],this_rrange,this_crange,false);
                    twist_storage(i,j) = this_twist;
                    rrange_storage(i,j,:) = permute(this_rrange,[1,3,2]);
                    crange_storage(i,j,:) = permute(this_crange,[1,3,2]);
                end
            end
            
            rrangeax = rranges(:,1);
            crangeax = cranges(:,1);
            
            if make_plots
                figure;
                imagesc(crangeax,rrangeax,twist_storage);
                title('Twist angles');
                axis equal
                ylabel('rrange lower limit')
                xlabel('crange lower limit')
                cbh = colorbar;
                ylabel(cbh,'Local moire twist angle (degrees)');
                set(gca,'yDir','normal');
            end
        end
        
        
        
        
        
        %%% Triangulation methods (Maddie, you may find these useful)
        % NPK revised 05/24/2020 to allow manual triangulation for dataset
        % 26, which does not triangulate correctly owing to distortion.
        function triangulateAA(obj,AA_type,use_manual_triangulation,crop_region)
            if nargin < 4
                crop_region = [];
            end
            if nargin < 3
                use_manual_triangulation = false;
            end
            if AA_type == 1  % circlefit
                warning('Reforming of axes could cause problems here');
                AA_cents = obj.AA_circlefit_values(:,1:2);
            elseif AA_type == 2  % Gaussian ellipse
                AA_cents = obj.AAgaussianfit_FWHM_ellipses(:,4:5);
            elseif AA_type == 3
                AA_cents = obj.AA_Gaussian_Circle_Fit_Params(:,1:2);
            else  % untrimmed AA centers.
                AA_cents = obj.untrimmed_AA_centers;
            end
            
            if ~isempty(crop_region)
                AA_cents(AA_cents(:,1)<crop_region{1}(1) | AA_cents(:,1)>crop_region{1}(2),:) = [];
                AA_cents(AA_cents(:,2)<crop_region{2}(1) | AA_cents(:,2)>crop_region{2}(2),:) = [];
            end
            
            if use_manual_triangulation
                disp('Beginning manual definition of triangulation connectivity.');
                % Recall that the following operations are in nm. Get
                % back to pixels for the triangulation itself.
                [f,ax] = obj.makeCustomDisplacementColorPlot([],[],[],[],1,0,0,0);
                hold on
                AA_cents_toplot = AA_cents*obj.scan_stepsize;
                scatter(AA_cents_toplot(:,1),AA_cents_toplot(:,2),100,'filled','g');
                count = 1;
                connectivity = zeros(0,3);
                while true
                    fprintf('Please click on three vertices for the %dth triangle. Clicks will snap to vertices.\n',count);
                    triangle_coords = ginput(3);
                    triangle_xcoords = triangle_coords(:,1)./obj.scan_stepsize;
                    triangle_ycoords = triangle_coords(:,2)./obj.scan_stepsize;
                    vertex_xcoords = AA_cents(:,1)';
                    vertex_ycoords = AA_cents(:,2)';
                    xdists = triangle_xcoords - vertex_xcoords;
                    ydists = triangle_ycoords - vertex_ycoords;
                    dists = sqrt(xdists.^2 + ydists.^2);
                    [~,vertex_indices] = min(dists,[],2);
                    connectivity(end+1,:) = vertex_indices;
                    count = count + 1;
                    yn = input('Do you wish to add another triangle? 1/0');
                    if ~yn
                        break
                    end
                end
                % This is actually not a delaunay triangulation, but keep
                % the DT nomenclature. The associated triangulation class
                % appears to possess all of the relevant methods, just with
                % manual definition of the adjacencies.
                DT = triangulation(connectivity,AA_cents);
            else
                DT = delaunayTriangulation(AA_cents);
            end
            obj.AAtriangulation = DT;
        end
        
        
        
        
        
        
        % Cleaning up the automatic triangulation, for samples such as the
        % tear, where there are too many triangles to be registered by
        % hand.
        function repairTriangulation(obj,ask)
            % for the triangulation error simulations
            if nargin < 2
                ask = true;
            end
            if isempty(ask)
                ask = true;
            end
            if isempty(obj.DSC_fit_storage)
                superimpose = false;
            else
                superimpose = true;
            end
            DT = obj.AAtriangulation;
            v = DT.Points;
            con = DT.ConnectivityList;
            % Repair on the basis of the tear mask, if it exists.
            if ~isempty(obj.tear_mask)
                tear_edge_mask = boundarymask(obj.tear_mask) & obj.tear_mask;
                [rvals,cvals] = find(tear_edge_mask);
                P = [cvals,rvals];
                P_toplot = (P-1)*obj.scan_stepsize;
                obj.plotFacetedTwistAngles([0,inf],parula,superimpose);
                hold on;
                scatter(P_toplot(:,1),P_toplot(:,2),'filled');
                title('Initial triangulation');
                ID = DT.pointLocation(P);
                ID(isnan(ID)) = [];
                uniqueIDs = unique(ID);
                % Because we cannot assign to the delaunayTriangulation, we 
                % will build a new triangulation (not a delaunay one) that
                % has removed the desired ids.
                con(uniqueIDs,:) = [];
                
            end
            
            % Repair by discarding all triangles with obtuse angles
            deleteinds = [];
            for i = 1:size(con,1)
                these_inds = con(i,:);
%                 these_coords = v(these_inds,:);
                v1 = v(these_inds(1),:);
                v2 = v(these_inds(2),:);
                v3 = v(these_inds(3),:);
                % angle between [v1v2 and v1v3]
                vec1 = v2-v1;
                vec2 = v3-v1;
                angletest1 = dot(vec1,vec2);
                % angle between [v2v1 and v2v3]
                vec1 = v1-v2;
                vec2 = v3-v2;
                angletest2 = dot(vec1,vec2);
                % angle between [v3v1 and v3v2]
                vec1 = v1-v3;
                vec2 = v2-v3;
                angletest3 = dot(vec1,vec2);
                test = angletest1 < 0 || angletest2 < 0 || angletest3 < 0;
                if test 
                    deleteinds(end+1) = i;
                end
            end
            con(deleteinds,:) = [];
            
            T = triangulation(con,v);
            obj.AAtriangulation = T;
            obj.setFacetTwistAngles();  % because the list has changed.
            
            obj.plotFacetedTwistAngles([0,inf],parula,superimpose);
            title('Triangulation after automated repair');
            
            % Now ask the user for input.
            f = gcf;
            if ask
                tf = input('Do you wish to perform manual repair? 1/0');
            else
                tf = 0;
            end
            if tf
                %                 disp('Click inside any triangulation facets to delete. Clicks will snap to facets.');
                deleteinds = zeros(0,1);
                figure(f);
                tf2 = input('Do you with to delete facets?');
                if tf2
                    while true
                        fprintf('Please click in any facet to remove it from the triangulation.\n');
                        click_coords = ginput(1);
                        click_xcoords = click_coords(:,1)./obj.scan_stepsize + 1;
                        click_ycoords = click_coords(:,2)./obj.scan_stepsize + 1;
                        ID = T.pointLocation([click_xcoords,click_ycoords]);
                        deleteinds(end+1,:) = ID;
                        yn = input('Do you wish to delete another triangle? 1/0');
                        if ~yn
                            break
                        end
                    end
                    if ~isempty(deleteinds)
                        con(deleteinds,:) = [];
                        T = triangulation(con,v);
                        obj.AAtriangulation = T;
                        obj.setFacetTwistAngles();  % because the list has changed.
                        obj.plotFacetedTwistAngles([0,inf],parula,superimpose);
                        title('Triangulation after manual facet deletion');
                    end
                end
                
                % Perform manual vertex addition
                tf3 = input('Do you wish to add facets?');
                if tf3
                    hold on;
                    v_toplot = (v-1)*obj.scan_stepsize;
                    scatter(v_toplot(:,1),v_toplot(:,2),100,'filled','g');
                    while true
                        fprintf('Please click on three vertices to form a triangle. Clicks will snap to vertices.\n');
                        click_coords = ginput(3);
                        triangle_xcoords = click_coords(:,1)./obj.scan_stepsize + 1;
                        triangle_ycoords = click_coords(:,2)./obj.scan_stepsize + 1;
                        vertex_xcoords = v(:,1)';
                        vertex_ycoords = v(:,2)';
                        xdists = triangle_xcoords - vertex_xcoords;
                        ydists = triangle_ycoords - vertex_ycoords;
                        dists = sqrt(xdists.^2 + ydists.^2);
                        [~,vertex_indices] = min(dists,[],2);
                        con(end+1,:) = vertex_indices;
                        yn = input('Do you wish to add another triangle? 1/0');
                        if ~yn
                            break
                        end
                    end
                end
                T = triangulation(con,v);
                obj.AAtriangulation = T;
                obj.setFacetTwistAngles();  % because the list has changed.
                obj.plotFacetedTwistAngles([0,inf],parula,true);
                title('Triangulation after manual facet repair');
            end
        end
        
        
        
        
        
        
        % Computes to the ith vertex in the untrimmed AA vertex list at
        % present. May be able to do this with trimmed later. Attached
        % vertices and attached triangles both come out as indices, to be
        % referenced against the list for triangulation computation.
        % (Untrimmed AA centers, currently)
        %
        % See code at the following link, adapted here:
        % https://www.mathworks.com/matlabcentral/answers/163115-how-do-i-get-the-neighbors-of-a-vertex-in-a-delaunay-triangulation
        function [attached_vertices,attached_triangles] = getNearestVertices(obj,idx)
            x = obj.untrimmed_AA_centers(:,1);
            y = obj.untrimmed_AA_centers(:,2);
            DT = obj.AAtriangulation;
            attached_triangles = DT.vertexAttachments(idx);
%             num_triangles_attached = numel(attached_triangles);
%             for i = 1:num_triangles_attached
                vertices_i = DT.ConnectivityList(attached_triangles{1},:);
                attached_vertices = setdiff(unique(vertices_i), idx);
%             end
            
            % Lets Check this
%             for j = 1: 5: size(x,1)
                figure
                triplot(DT);
                hold on;
                plot(x(idx),y(idx),'rx')
                plot(x(attached_vertices),y(attached_vertices),'ro')
                hold off;
        end
        
        
        
        
        
        function setFacetTwistAngles(obj)
            obj.reformAxes();
            DT = obj.AAtriangulation;
            con = DT.ConnectivityList;
            v = DT.Points;
            ntriangles = size(con,1);
            obj.triangle_sidelengths = zeros(ntriangles,3);
            obj.triangle_moire_angle = zeros(ntriangles,1);
            
            [ ~, ~, hexagon_lattice_constant ] = getDSCBasisVectors();
            a = hexagon_lattice_constant;
            for i = 1:ntriangles
                these_vertices = con(i,:);
                v1 = v(these_vertices(1),:);
                v2 = v(these_vertices(2),:);
                v3 = v(these_vertices(3),:);
                % now that we have the coordinates, get all side lengths.
                l1 = sum((v1-v2).^2).^0.5*obj.scan_stepsize;
                l2 = sum((v1-v3).^2).^0.5*obj.scan_stepsize;
                l3 = sum((v2-v3).^2).^0.5*obj.scan_stepsize;
                obj.triangle_sidelengths(i,:) = [l1,l2,l3];
                meanl = mean([l1,l2,l3]);
                % NPK note 05/27/2020: the first formula I used below
                % implicitly employed the small angle approximation. The
                % following formula removes that approximation. However,
                % the small angle approximation is better than one part in
                % ten thousand for our datasets, so it gives the same
                % numbers as before.
%                 obj.triangle_moire_angle(i) = rad2deg((a/10)/meanl);
                ause = a/10;
                radangle = 2*asin(ause/(2*meanl));
                obj.triangle_moire_angle(i) = rad2deg(radangle);
            end
            
            % Triangles on the boundary are suspect. Trim these out.
            
            %edges 
            
        end
        
        
        
        
        
        
        % Requires the above methods to be run
        % Excludes all triangulation facets that do not fall within a
        % particular angle range, which is the same as how we make the
        % plots above. Use plotting to figure out what this value is.
        %
        % Uses the method in the SI of "Maximized electron interactions at
        % the magic angle in twisted bilayer graphene", Kerelsky et al.,
        % Nature 2019.
        %
        % The return variable nconvergedstor is a check to make sure the
        % fitting function is optimizing correctly. 
        function [average_rms_residual,nconvergedstor] = fitUniaxialStrain(obj,average_angle_range,angle_guess)
            USE_OLD_FUNCTION = false;
            if isempty(average_angle_range)
                average_angle_range = [0,inf]; 
            end
            if nargin < 3
                angle_guess = 1;
            end
            triangle_lengths = obj.triangle_sidelengths;
            angle_values = obj.triangle_moire_angle;
            triangle_lengths(angle_values < average_angle_range(1) | angle_values > average_angle_range(2),:) = nan;
            n_facets = size(triangle_lengths,1);
            n_triangles = n_facets - nnz(isnan(triangle_lengths(:,1)));
            uniaxial_moire_angle_storage = zeros(n_facets,1);
            uniaxial_strain_angle_storage = zeros(n_facets,1);
            uniaxial_strain_percent_storage = zeros(n_facets,1);
            options.Display = 'off';
            nconvergedstor = zeros(n_facets,1);
            bestrmsstor = zeros(n_facets,1);
            bestresidstor = zeros(n_facets,3);
            for i = 1:n_facets
                if isnan(triangle_lengths(i,1))
                    uniaxial_moire_angle_storage(i) = nan;
                    uniaxial_strain_angle_storage(i) = nan;
                    uniaxial_strain_percent_storage(i) = nan;
                    continue
                end
                if mod(i,10) == 0
                    fprintf('Fitting Moire triangle %d of %d.\n',i,n_facets);
                end
%             [ moire_wavelengths ] = uniaxialStrainPredFun( theta_moire, theta_strain, epsilon_strain_percent );
                
                if USE_OLD_FUNCTION
                    these_expt_wavelengths = sort(triangle_lengths(i,:));
                    uniaxialStrainResidFun = @(p) these_expt_wavelengths - uniaxialStrainPredFun(p(1),p(2),p(3));
                else
                    these_expt_wavelengths = triangle_lengths(i,:);
                    uniaxialStrainResidFun = @(p) these_expt_wavelengths - uniaxialStrainPredFun2(p(1),p(2),p(3));
                end
%                 init_guesses = [1,0,0.1;
% %                                 1,45,0.1;
%                                 1,90,0.1;
% %                                 1,135,0.1;
%                                 1,180,0.1;
% %                                 1,225,0.1;
%                                 1,270,0.1];
% %                                 1,315,0.1];
                init_guesses = [angle_guess,2,0.1;
                                angle_guess,45,0.1;
                                angle_guess,90,0.1;
                                angle_guess,135,0.1;
                                angle_guess,178,0.1];
%                                 angle_guess,225,0.1;
%                                 angle_guess,270,0.1];
%                                 angle_guess,315,0.1];
                %   Preliminary tests suggest complete convergence, so cut
                %   down on number of multistarts.
                rmsstor = zeros(size(init_guesses,1),1);
                inner_resid_stor = zeros(size(init_guesses,1),3);
                fitvalstor = zeros(size(init_guesses,1),3);
                if USE_OLD_FUNCTION
                    lb = [0,-360,0];
                    ub = [360,360,inf];
                else
                    lb = [0,0,0];
                    ub = [360,180,inf];
                end
                for j = 1:size(init_guesses,1)
                    this_inital_guess = init_guesses(j,:);
                    fitvals = lsqnonlin(uniaxialStrainResidFun,this_inital_guess,lb,ub,options);
                    rmsval = rms(uniaxialStrainResidFun(fitvals));
                    rmsstor(j) = rmsval;
                    residval = uniaxialStrainResidFun(fitvals);
                    inner_resid_stor(j,:) = residval;
                    fitvalstor(j,:) = fitvals;
                end
                [val,idx] = min(rmsstor);
                best_params = fitvalstor(idx,:);
                tol = 1e-4;  % criterion for it actually being the same convergence location.
                nconverged = nnz(abs(rmsstor - val) < tol);
                
                uniaxial_moire_angle_storage(i) = best_params(1);
                uniaxial_strain_angle_storage(i) = best_params(2);
                uniaxial_strain_percent_storage(i) = best_params(3);
                nconvergedstor(i) = nconverged;
                bestrms = rmsstor(idx);
                bestrmsstor(i) = bestrms;
                bestresidstor(i,:) = inner_resid_stor(idx,:);
            end
            obj.uniaxial_moire_angles = uniaxial_moire_angle_storage;
            obj.uniaxial_strain_angles = uniaxial_strain_angle_storage;
            obj.uniaxial_strain_percents = uniaxial_strain_percent_storage;
            obj.uniaxial_strain_residuals = bestresidstor;
            average_rms_residual = mean(bestrmsstor,'omitnan');
        end
        
        
        
        
        
        function [av_moire_angle,std_moire_angle,min_moire_angle,max_moire_angle,...
                av_strain_percent,std_strain_percent,min_strain_percent,max_strain_percent,copypaste] = ...
                getUniaxialStrainStatistics(obj,make_plots,moire_binedges,strainangle_binedges,strainval_binedges)
            if (nargin < 3), moire_binedges = []; end
            if (nargin < 4), strainangle_binedges = []; end
            if (nargin < 5), strainval_binedges = []; end
            av_moire_angle = mean(obj.uniaxial_moire_angles,'omitnan');
            std_moire_angle = std(obj.uniaxial_moire_angles,'omitnan');
            min_moire_angle = min(obj.uniaxial_moire_angles,[],'omitnan');
            max_moire_angle = max(obj.uniaxial_moire_angles,[],'omitnan');
            av_strain_percent = mean(obj.uniaxial_strain_percents,'omitnan');
            std_strain_percent = std(obj.uniaxial_strain_percents,'omitnan');
            min_strain_percent = min(obj.uniaxial_strain_percents,[],'omitnan');
            max_strain_percent = max(obj.uniaxial_strain_percents,[],'omitnan');
            copypaste = [av_moire_angle;std_moire_angle;min_moire_angle;max_moire_angle;av_strain_percent;std_strain_percent;min_strain_percent;max_strain_percent];
            
            if make_plots
                use_eps = obj.uniaxial_strain_percents;
                use_eps(isnan(use_eps)) = [];
                figure; 
                if isempty(strainval_binedges)
                    histogram(use_eps);
                else
                    histogram(use_eps,strainval_binedges);
                end
                xlabel('Uniaxial strain %');
                ylabel('Counts');
                
                use_mangles = obj.uniaxial_moire_angles;
                use_mangles(isnan(use_mangles)) = [];
                figure; 
                if isempty(moire_binedges)
                    histogram(use_mangles);
                else
                    histogram(use_mangles,moire_binedges);
                end
                xlabel('Uniaxial-calculated Moire twist angle (degrees)');
                ylabel('Counts');
                
                use_sangles = obj.uniaxial_strain_angles;
                use_sangles(isnan(use_sangles)) = [];
                figure; 
                if isempty(strainangle_binedges)
                    histogram(use_sangles);
                else
                    histogram(use_sangles,strainangle_binedges);
                end
                xlabel('Uniaxial-calculated strain angle (degrees)');
                ylabel('Counts');
                
                % Correlation between strain value and strain angle
                figure; 
                scatter(use_sangles,use_eps);
                xlabel('Uniaxial strain angle (degrees)');
                ylabel('Uniaxial strain value (degrees)');
            end
        end
        
        
        
        
        
        % No need for angle_range because those rows of the data have
        % already been set to nans in the function computing the uniaxial
        % strain.
        function plotFacetedUniaxialStrain(obj,cmap,removeLines)
            strain_percents = obj.uniaxial_strain_percents;
            strain_angles = obj.uniaxial_strain_angles;
            DT = obj.AAtriangulation;
            con = DT.ConnectivityList;
            % Added 06/02/2020 so that we can do a superimposition
            p = DT.Points;
            p = (p-1)*obj.scan_stepsize;
            DTplot = triangulation(con,p);
            con = DTplot.ConnectivityList;
            v = DTplot.Points;
            ntriangles = size(con,1);
            if nargin < 2
                cmap = parula;
            end
            
            for q = 1:2
                if q == 1
                    mat = strain_percents;
                else
                    mat = strain_angles;
                end
                max_val = max(mat);
                min_val = min(mat);
                interpanchors = linspace(min_val,max_val,size(cmap,1));
                figure;
                if ~(removeLines)
                    triplot(DTplot);
                end
                hold on
                for i = 1:ntriangles
                    % get vertices for patch, figure out what color to use.
                    thisval = mat(i);
                    if isnan(thisval)
                        continue
                    end
                    c = zeros(1,3);
                    c(1) = interp1(interpanchors,cmap(:,1),thisval);
                    c(2) = interp1(interpanchors,cmap(:,2),thisval);
                    c(3) = interp1(interpanchors,cmap(:,3),thisval);
                    these_vertices = con(i,:);
                    v1 = v(these_vertices(1),:);
                    v2 = v(these_vertices(2),:);
                    v3 = v(these_vertices(3),:);
                    coords = vertcat(v1,v2,v3);
                    h = patch(coords(:,1),coords(:,2),c);
                    if removeLines
                        h.EdgeColor = 'none';
                    end
                end
                
                
                if removeLines
                    DTplot2 = [];
                    if ~isempty(obj.tear_mask)
                        % Separate triangulation into two components
                        pointsall = DTplot.Points;
                        [maskr,maskc] = find(obj.tear_mask);
                        leftpoint_rinds = pointsall(:,1)./obj.scan_stepsize < mean(maskc);
                        leftpoints = pointsall(leftpoint_rinds,:);
                        rightpoint_rinds = pointsall(:,1)./obj.scan_stepsize > mean(maskc);
                        rightpoints = pointsall(rightpoint_rinds,:);
                        leftedges = zeros(0,3);
                        rightedges = zeros(0,3);
                        leftpoint_rinds_nums = find(leftpoint_rinds);
                        rightpoint_rinds_nums = find(rightpoint_rinds);
                        for k = 1:size(DTplot.ConnectivityList,1)
                            thiscon = DTplot.ConnectivityList(k,:);
                            if nnz(nnz(thiscon == leftpoint_rinds_nums))
                                leftedges(end+1,:) = thiscon;
                            else
                                rightedges(end+1,:) = thiscon;
                            end
                        end
                        % Quite hacky but works
                        DTplot = triangulation(rightedges-size(leftpoints,1),rightpoints);
                        try
                            DTplot2 = triangulation(leftedges,leftpoints);
                        catch
                        end
                    end
                    F = freeBoundary(DTplot);
                    xpoints = DTplot.Points(:,1);
                    ypoints = DTplot.Points(:,2);
                    xplot = xpoints(F);
                    yplot = ypoints(F);
                    hold on;
                    plot(xplot,yplot,'Color',[0.8,0.8,0.8]);
                    if false && ~isempty(DTplot2)
                        F = freeBoundary(DTplot2);
                        xpoints = DTplot2.Points(:,1);
                        ypoints = DTplot2.Points(:,2);
                        xplot = xpoints(F);
                        yplot = ypoints(F);
                        hold on;
                        plot(xplot,yplot,'Color',[0.8,0.8,0.8]);
                    end
                end
                colormap(cmap);
                cbh = colorbar;
                caxis([min_val,max_val]);
                if q == 1
                    ylabel(cbh,'Calculated uniaxial strain (%)');
                    title('Uniaxial strain %');
                else
                    ylabel(cbh,'Calculated uniaxial strain angle (deg)');
                    title('Uniaxial strain angle');
                end
                axis equal
                set(gca,'yDir','normal');
            end
            
        end
        
        
        
        
        
        function plotFacetedTwistAngleHistogram(obj,moire_binedges)
            if nargin < 2
                moire_binedges = [];
            end
            use_mangles = obj.triangle_moire_angle;
            use_mangles(isnan(use_mangles)) = [];
            figure;
            if isempty(moire_binedges)
                histogram(use_mangles);
            else
                histogram(use_mangles,moire_binedges);
            end
            xlabel('Average sidelength-calculated Moire twist angle (degrees)');
            ylabel('Counts');
        end
        
        
        
        
        
        
        function [average_angle,std_angle,max_angle,min_angle,n_triangles_used,angles_used] = plotFacetedTwistAngles(obj,anglerange,cmap,superimpose,use_uniaxial_angles,removelines,caxislimits)
            if isempty(anglerange)
                anglerange = [0,inf];
            end
            if nargin < 4
                superimpose = true;
            end
            if nargin < 5
                use_uniaxial_angles = false;
            end
            if nargin < 6
                removelines = false;
            end
            if use_uniaxial_angles
                moirean_touse = obj.uniaxial_moire_angles;
            else
                moirean_touse = obj.triangle_moire_angle;
            end
            if nargin > 1
                moirean_touse(moirean_touse < anglerange(1)) = nan;
                moirean_touse(moirean_touse > anglerange(2)) = nan;
            end
            DT = obj.AAtriangulation;
            con = DT.ConnectivityList;
            % Added 06/02/2020 so that we can do a superimposition
            p = DT.Points;
            p = (p-1)*obj.scan_stepsize;
            DTplot = triangulation(con,p);
            if superimpose 
                [figh,axh] = obj.makeCustomDisplacementColorPlot();
                figure(figh);
                hold on
            else
                figure;
            end
            v = DTplot.Points;
            ntriangles = size(con,1);
            if nargin < 3
                cmap = parula;
            end
            if nargin < 7
                caxislimits = [];
            end
            if isempty(caxislimits)
                max_angle = max(moirean_touse);
                min_angle = min(moirean_touse);
            else
                max_angle = max(caxislimits);
                min_angle = min(caxislimits);
            end
            interpanchors = linspace(min_angle,max_angle,size(cmap,1));
            if ~(removelines)
                triplot(DTplot);
            end
            hold on
            for i = 1:ntriangles
                % get vertices for patch, figure out what color to use.
                thisangle = moirean_touse(i);
                if isnan(thisangle)
                    continue
                end
                c = zeros(1,3);
                c(1) = interp1(interpanchors,cmap(:,1),thisangle);
                c(2) = interp1(interpanchors,cmap(:,2),thisangle);
                c(3) = interp1(interpanchors,cmap(:,3),thisangle);
                these_vertices = con(i,:);
                v1 = v(these_vertices(1),:);
                v2 = v(these_vertices(2),:);
                v3 = v(these_vertices(3),:);
                coords = vertcat(v1,v2,v3);
                h = patch(coords(:,1),coords(:,2),c);
                if removelines
                    h.EdgeColor = 'none';
                end
                if superimpose
                    h.FaceAlpha = 0.9;
                end
            end
            if removelines
                DTplot2 = [];
                if ~isempty(obj.tear_mask)
                    % Separate triangulation into two components
                    pointsall = DTplot.Points;
                    [maskr,maskc] = find(obj.tear_mask);
                    leftpoint_rinds = pointsall(:,1)./obj.scan_stepsize < mean(maskc);
                    leftpoints = pointsall(leftpoint_rinds,:);
                    rightpoint_rinds = pointsall(:,1)./obj.scan_stepsize > mean(maskc);
                    rightpoints = pointsall(rightpoint_rinds,:);
                    leftedges = zeros(0,3);
                    rightedges = zeros(0,3);
                    leftpoint_rinds_nums = find(leftpoint_rinds);
                    rightpoint_rinds_nums = find(rightpoint_rinds);
                    for q = 1:size(DTplot.ConnectivityList,1)
                        thiscon = DTplot.ConnectivityList(q,:);
                        if nnz(nnz(thiscon == leftpoint_rinds_nums))
                            leftedges(end+1,:) = thiscon;
                        else
                            rightedges(end+1,:) = thiscon;
                        end
                    end
                    % Quite hacky but works
                    DTplot = triangulation(rightedges-size(leftpoints,1),rightpoints);
                    try
                        DTplot2 = triangulation(leftedges,leftpoints);
                    catch
                    end
                end
                F = freeBoundary(DTplot); 
                xpoints = DTplot.Points(:,1);
                ypoints = DTplot.Points(:,2);
                xplot = xpoints(F);
                yplot = ypoints(F);
                hold on;
                plot(xplot,yplot,'Color',[0.8,0.8,0.8]);
                if false && ~isempty(DTplot2)
                    F = freeBoundary(DTplot2);
                    xpoints = DTplot2.Points(:,1);
                    ypoints = DTplot2.Points(:,2);
                    xplot = xpoints(F);
                    yplot = ypoints(F);
                    hold on;
                    plot(xplot,yplot,'Color',[0.8,0.8,0.8]);
                end
            end
            cbh = colorbar;
            colormap(cmap);
            caxis([min_angle,max_angle]);
            ylabel(cbh,'Calculated Moire twist angle (deg)');
            axis equal
            set(gca,'yDir','normal');
            title('Geometric twist angle analysis');
            xlabel('x (nm)');
            ylabel('y (nm)');
            
            moirean_touse(isnan(moirean_touse)) = [];
            average_angle = mean(moirean_touse);
            std_angle = std(moirean_touse);
            n_triangles_used = numel(moirean_touse);
            angles_used = moirean_touse;
        end
            
            
         
        
        
        % This function sets the instance variable
        % obj.moire_angle_estimates, which contains an estimate for the
        % effective moire twist angle at each pixel in the real space
        % image. For large twist angles, this is best produced by the
        % moving average method given above. For small twist angles, this
        % best produced by the triangulation method. In both cases, the 
        % twist angle estimate will be blurred by a filter so that sharp
        % edges from the estimation method do not adversely impact the
        % computed strain fields. The unblurred result will be called 
        % obj.moire_angle_estimates_raw.
        %
        % Nathanael Kazmierczak, 05/22/2020
        function computeMoireAngleEstimates(obj,method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input)
            if strcmp(method,'outside input')
                angles = angle_input*ones(obj.datacube_size(3:4));
                obj.moire_angle_estimates_raw = angles;
                obj.moire_angle_estimates = angles;
            elseif strcmp(method,'single window')
                crop_margin = crop_val_pixels;
                rrange = [];
                crange = [];
                make_plots = true;
                twist_angle = obj.fitMoireAngle(crop_margin,rrange,crange,make_plots);
                angles = twist_angle*ones(obj.datacube_size(3:4));
                obj.moire_angle_estimates_raw = angles;
                obj.moire_angle_estimates = angles;
            elseif strcmp(method,'moving window')
                rvals = numel(obj.xaxis);
                cvals = numel(obj.yaxis);
                rrange_global = [crop_val_pixels,rvals-crop_val_pixels];
                crange_global = [crop_val_pixels,cvals-crop_val_pixels];
                [twist_storage,rrange_storage,crange_storage] = obj.fitMovingMoireAngle(rrange_global,crange_global,window_size,window_increment);
                % Each pixel will now get assigned to the window it is
                % closest to.
                centerr = mean(rrange_storage,3);
                centerc = mean(crange_storage,3);
                rcents = centerr(:)';
                ccents = centerc(:)';
                
                rbase = 1:obj.datacube_size(3);
                cbase = 1:obj.datacube_size(4);
                [cspace,rspace] = meshgrid(rbase,cbase);
                
                rdiffs = rspace(:) - rcents;
                cdiffs = cspace(:) - ccents;
                dists = sqrt(rdiffs.^2 + cdiffs.^2);
                [~,idxs] = min(dists,[],2);
                angle_vals_lin = twist_storage(idxs);
                angle_vals = reshape(angle_vals_lin,size(rspace));
                figure; imagesc(angle_vals); colorbar; set(gca,'yDir','normal');
                obj.moire_angle_estimates_raw = angle_vals;
                stringcell = {'moving average circle','xdisp',15};
                filterstruct = buildFilterStruct(stringcell);
                [ angles_filt,~ ] = filterDisplacement( angle_vals,[],filterstruct,0,[] );
                figure; imagesc(angles_filt); colorbar; set(gca,'yDir','normal');
                obj.moire_angle_estimates = angles_filt;
            elseif strcmp(method,'triangulation')
                if strcmp(AA_type_string,'AA Gaussian ellipse')
                    AA_type = 2;
                elseif strcmp(AA_type_string,'AA mask circle')
                    AA_type = 1;
                elseif strcmp(AA_type_string,'AA untrimmed')
                    AA_type = 0;
                else 
                    error('Please choose a valid option for AA type in computing the triangulated moire estimate.');
                end
                obj.triangulateAA(AA_type);
                obj.setFacetTwistAngles();
            else
                error('Please select either ''moving window'' or ''triangulation'' for the moire angle estimate method!');
            end
            
        end
        
        
        
        
        
        
        % Strain mapping mask functions
        % Buffer should be returned as full matrix so that it can be
        % trimmed to fit the reliable strain mapping zone. Be sure to plot
        % mask.
        % Use range within ring as a reliability metric to figure out when
        % we have left the AA region.
        % Implement a circular buffer first, but an elliptical buffer might
        % be good later.
        % Both statistics outputs are for the smallest buffer size, at the
        % moment.
        function [averages_mat,stds_mat,n_values,names,stderrs,nconn] = getStrainInAARingedBufferMasks(obj,...
                AA_center_type,buffer_type,radius_values,saveplotflag,AA_crop_region,strain_crop_region)
            if nargin < 6
                AA_crop_region = [];
                strain_crop_region = [];
            end
            
            switch AA_center_type
                case 'untrimmed AA'
                    AA_cents = obj.untrimmed_AA_centers;
                case 'Mask circle AA fit'
                    AA_cents = obj.AA_circlefit_values(:,1:2);
                case 'Gaussian ellipse AA fit'
                    AA_cents = obj.AAgaussianfit_FWHM_ellipses(:,4:5);
                otherwise 
                    error('Please choose a valid AA_center_type');
            end
            
            if ~isempty(AA_crop_region)
                AA_cents(AA_cents(:,1)<AA_crop_region{1}(1) | AA_cents(:,1)>AA_crop_region{1}(2),:) = [];
                AA_cents(AA_cents(:,2)<AA_crop_region{2}(1) | AA_cents(:,2)>AA_crop_region{2}(2),:) = [];
            end
            
            if strcmp(buffer_type,'circular')
                xbase = obj.xaxis;
                ybase = obj.yaxis;
                [xspace,yspace] = meshgrid(xbase,ybase);
                [rs,cs] = size(xspace);
                [tr,tc] = size(obj.exx);
                crop_value = (rs-tr)/2;
                nr = numel(radius_values);
                radius_inclusive_masks = false(rs,cs,nr);
                centers = (AA_cents-1)*obj.scan_stepsize;  % because all of the AA centers are originally in pixels
                ncent = size(centers,1);
                % The -1 acccounts that the axes are zero indexed while
                % the pixels are 1 indexed.
                for i = 1:nr
                    this_radius = radius_values(i);
                    tfstor = false(rs,cs);
                    for j = 1:ncent
                        this_cent = centers(j,:);
                        [tf] = isInCircle(xspace,yspace,this_cent(1),this_cent(2),this_radius); 
                        tfstor = tf | tfstor;
                    end
                    radius_inclusive_masks(:,:,i) = tfstor;
                end
                
                % Convert the inclusive masks to a layered ringed buffer
                radius_exclusive_masks = false(size(radius_inclusive_masks));
                TESTPLOTS = false;
                for k = 1:nr
                    % Make sure summation through a logical works
                    has_been_allocated = sum(radius_exclusive_masks,3);
                    this_mask = radius_inclusive_masks(:,:,k) & ~has_been_allocated;
                    radius_exclusive_masks(:,:,k) = this_mask;
                    
                    if TESTPLOTS
                        figure
                        subplot(1,2,1);
                        imagesc(obj.xaxis,obj.yaxis,radius_inclusive_masks(:,:,k)); axis equal; set(gca,'yDir','normal');
                        subplot(1,2,2);
                        imagesc(obj.xaxis,obj.yaxis,radius_exclusive_masks(:,:,k)); axis equal; set(gca,'yDir','normal');
                    end
                end
                
                % Trim the masks
                % Added 11/12/2020: This is probably the best way to do it
                % generally, but leave what was previously in place for
                % consistency
                if ~isempty(strain_crop_region)
                    crop_value = obj.trim_value;
                end
                exclusive_masks_trimmed = false([size(trimArray(radius_exclusive_masks(:,:,1),crop_value)),nr]);
                inclusive_masks_trimmed = false([size(trimArray(radius_inclusive_masks(:,:,1),crop_value)),nr]);
                for k = 1:nr
                    exclusive_masks_trimmed(:,:,k) = trimArray(radius_exclusive_masks(:,:,k),crop_value);
                    inclusive_masks_trimmed(:,:,k) = trimArray(radius_inclusive_masks(:,:,k),crop_value);
                end
                
                % Added 11/12/2020: 
                if ~isempty(strain_crop_region)
                    exclusive_masks_trimmed = false(numel(strain_crop_region{1}),numel(strain_crop_region{2}),nr);
                    inclusive_masks_trimmed = false(numel(strain_crop_region{1}),numel(strain_crop_region{2}),nr);
%                     inclusive_masks_trimmed = false([size(trimArray(radius_inclusive_masks(:,:,1),crop_value)),nr]);
                for k = 1:nr
                    exclusive_masks_trimmed_2(:,:,k) = trimArray(radius_exclusive_masks(:,:,k),crop_value);
                    inclusive_masks_trimmed_2(:,:,k) = trimArray(radius_inclusive_masks(:,:,k),crop_value);
                    exclusive_masks_trimmed(:,:,k) = exclusive_masks_trimmed_2(strain_crop_region{1},strain_crop_region{2},k);
                    inclusive_masks_trimmed(:,:,k) = inclusive_masks_trimmed_2(strain_crop_region{1},strain_crop_region{2},k);
                end
                end
                
                % Apply the masks to the different components of the strain
                % field.
                matcell = {obj.exx,obj.exy,obj.eyx,obj.eyy,obj.gxy,obj.fixed_body_rotation,obj.principal_strain_1};
                legendcell = {'exx','exy','eyx','eyy','gxy','fixed_body_rotation','principal strain max'};
                nmats = numel(matcell);  % what we have stored and of interest at the moment
                storage_mat = zeros(nr,5,nmats);  % mean, median, std, max, min
                for i = 1:nr
                    for j = 1:nmats
                        this_mat = matcell{j};
                        this_emask = exclusive_masks_trimmed(:,:,i);
                        these_vals = this_mat(this_emask);
                        storage_mat(i,:,j) = [mean(these_vals),median(these_vals),std(these_vals),max(these_vals),min(these_vals)];
                    end
                end
                
                % Make plots
                for q = 1:nmats
                    figh = figure
                    plot(radius_values,storage_mat(:,1,q),'-o');
                    hold on
                    plot(radius_values,storage_mat(:,4,q),'-o');
                    plot(radius_values,storage_mat(:,5,q),'-o');
                    legend('Mean','Max','Min');
                    xlabel('AA ringed buffer radius size (exclusive) (nm)');
                    ylabel('Strain %');
                    title(sprintf('%s, %s',obj.filename,legendcell{q}),'Interpreter','none');
                    
                    if saveplotflag
                        currentd = pwd;
                        cd(obj.saved_plots_folderpath);
                        if ~exist(obj.saved_plots_foldername,'dir')
                            mkdir(obj.saved_plots_foldername);
                        end
                        cd(obj.saved_plots_foldername);
                        saveas(figh,sprintf('AAeringedBuffer_%s.png',legendcell{q}));
                        savefig(figh,sprintf('AAeringedBuffer_%s',legendcell{q}));
                        cd(currentd);
                    end
                end
            end
            
            if nargout > 4  % need standard errors, so segment mask into connected components
               % This had better be a matrix: inclusive_masks_trimmed
                mask_segmentation = bwconncomp(inclusive_masks_trimmed);
                nconn = mask_segmentation.NumObjects;
                component_storage_mat = zeros(nmats,nconn);  % mean only, different conn comps in cols
                for q = 1:nconn
                    this_mask = false(size(inclusive_masks_trimmed));
                    this_mask(mask_segmentation.PixelIdxList{q}) = true;
                    for s = 1:nmats
                        % matcell has not changed, can run these
                        % calculations for each matrix separately.
                        this_mat = matcell{s};
                        these_vals = this_mat(this_mask);
                        component_storage_mat(s,q) = mean(these_vals);
                    end
                end
                % get std for each type of matrix (each row).
                stdevs_comp = std(component_storage_mat,0,2); % second dimension, default normalization
                stderrs = stdevs_comp / sqrt(nconn);
            end
            
            
            % Figure out the average of each type of strain in the center
            % of the AA domains. Use only the smallest exclusive mask.
            averages_mat = permute(storage_mat(1,[1,3],:),[3,2,1]);
            stds_mat = averages_mat(:,2);
            averages_mat = averages_mat(:,1);
            n_values = nnz(exclusive_masks_trimmed(:,:,1));
            names = legendcell';
        end
        
        
        
        
        % Helper function for the strain buffering methods
        function [untrimmed_mask,trimmed_mask] = getAACircleBufferMasks(obj,AA_center_type,radius_value)
            switch AA_center_type
                case 'untrimmed AA'
                    AA_cents = obj.untrimmed_AA_centers;
                case 'Mask circle AA fit'
                    AA_cents = obj.AA_circlefit_values(:,1:2);
                case 'Gaussian ellipse AA fit'
                    AA_cents = obj.AAgaussianfit_FWHM_ellipses(:,4:5);
                case 'Gaussian circle AA fit'
                    AA_cents = obj.AA_Gaussian_Circle_Fit_Params(:,1:2);
                otherwise
                    error('Please choose a valid AA_center_type');
            end
            
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            [xspace,yspace] = meshgrid(xbase,ybase);
            [rs,cs] = size(xspace);
            [tr,tc] = size(obj.exx);
            crop_value = (rs-tr)/2;
            centers = (AA_cents-1)*obj.scan_stepsize;  % because all of the AA centers are originally in pixels
            ncent = size(centers,1);
            % The -1 acccounts that the axes are zero indexed while
            % the pixels are 1 indexed.
            tfstor = false(rs,cs);
            for j = 1:ncent
                this_cent = centers(j,:);
                [tf] = isInCircle(xspace,yspace,this_cent(1),this_cent(2),radius_value);
                tfstor = tf | tfstor;
            end
            untrimmed_mask = tfstor;
            trimmed_mask = trimArray(tfstor,crop_value);
            
            TEST_FIGURES = true;
            if TEST_FIGURES
                figure;
                imagesc(obj.xaxis,obj.yaxis,untrimmed_mask);set(gca,'yDir','normal');axis equal;
                title('Untrimmed mask');
                figure
                imagesc(obj.xaxis,obj.yaxis,trimmed_mask);set(gca,'yDir','normal');axis equal;
                title('Trimmed mask');
            end
        end
        
        
        
        
        % Buffer SP regions using the BWdist function 
        % Recall that a SPdirection value of 1 is the purple saddle points,
        % which is how we set the x axis for the plotting by default.
        function [untrimmed_mask] = makeSaddlePointMask(obj,SPdirection,buffer_distance,useStrainFilterMatchedSP)
            if useStrainFilterMatchedSP
                storage_cell = obj.SP_registration_unique_strainfiltermatched;
            else
                storage_cell = obj.soliton_walls_unique;
            end
            if strcmp(SPdirection,'all')
                SP_mask = full(storage_cell{1}) | full(storage_cell{2}) | full(storage_cell{3});
            else
                SP_mask = full(storage_cell{SPdirection});
            end
            distance_mat = bwdist(SP_mask);
            addpoints = distance_mat <= buffer_distance;
            untrimmed_mask = SP_mask | addpoints;
            
            TEST_FIGURES = false;
            if TEST_FIGURES
                figure;
                imagesc(obj.xaxis,obj.yaxis,untrimmed_mask);set(gca,'yDir','normal');axis equal;
                title('Untrimmed mask');
            end
        end
        
        
        
        
        
        
        % A single shot function. We can also make a multi-ringed buffer
        % later. This function removes the AA regions.
        %
        % NPK modification on 06/17/2020: add bwconncomp to segment the
        % mask into all of the connected components. Then use the values in
        % the connected components to obtain the standard errors.
        %
        % nconn gives the number of connected components used for getting
        % the standard error.
        function [values,stds,npixels,names,stderrs,nconn] = getStrainInSPBufferMask(obj,SPdirection,SPbufferdist,AAremovalradius,AAcentertype,useStrainFilterMatchedSP,strain_crop_range)
            if nargin < 7
                strain_crop_range = [];
            end
            untrimmed_SP_masks = cell(0,1);
            if strcmp(SPdirection,'all')
                untrimmed_SP_masks{1} = obj.makeSaddlePointMask(1,SPbufferdist,useStrainFilterMatchedSP);
                untrimmed_SP_masks{2} = obj.makeSaddlePointMask(2,SPbufferdist,useStrainFilterMatchedSP);
                untrimmed_SP_masks{3} = obj.makeSaddlePointMask(3,SPbufferdist,useStrainFilterMatchedSP);
            else
                untrimmed_SP_masks{1} = obj.makeSaddlePointMask(SPdirection,SPbufferdist,useStrainFilterMatchedSP);
            end
            [untrimmed_AA_mask,trimmed_AA_mask] = obj.getAACircleBufferMasks(AAcentertype,AAremovalradius);
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            [xspace,yspace] = meshgrid(xbase,ybase);
            [rs,cs] = size(xspace);
            [tr,tc] = size(obj.exx);
            crop_value = (rs-tr)/2;
            SP_only_masks = cell(0,1);
            for i = 1:numel(untrimmed_SP_masks)  % Need to form all of the needed masks.
                this_SP_only_umask = untrimmed_SP_masks{i} & ~untrimmed_AA_mask;
                
                if ~isempty(strain_crop_range)
                    SP_only_masks{i} = trimArray(this_SP_only_umask,obj.trim_value);
                    SP_only_masks{i} = SP_only_masks{i}(strain_crop_range{1},strain_crop_range{2});
                else
                    SP_only_masks{i} = trimArray(this_SP_only_umask,crop_value);
                end
                TEST_FIGURES = true;
                if TEST_FIGURES
                    figure;
                    imagesc(obj.xaxis,obj.yaxis,SP_only_masks{i});set(gca,'yDir','normal');axis equal;
                    title('Trimmed SP mask');
                end
            end
            
            % Modification 06/23/2020: loop over all possible masks,
            % because we may need to perform tensor rotation between them.
            matcell_notrot = {obj.exx,obj.exy,obj.eyx,obj.eyy,obj.gxy,obj.fixed_body_rotation,obj.principal_strain_1,obj.principal_strain_diff};
            legendcell = {'exx','exy','eyx','eyy','gxy','fixed body rotation','principal strain max component','maximum (principal) shear'};
            nmats = numel(matcell_notrot);  % what we have stored and of interest at the moment
            storage_mat = zeros(nmats,5);  % mean, median, std, max, min
            
            % compute the value for each quantity
            if numel(SP_only_masks) > 1
                storage_cell = cell(nmats,numel(SP_only_masks));  % just storing values, not mean, median, etc
                for i = 1:numel(SP_only_masks)
                    %                     sind = (i-1)*nmats + 1;
                    %                     eind = i*nmats;
                    % Include tensor rotation of the strain fields for the
                    % second and third soliton walls.
                    if i == 1
                        % Assumption for the time being is that the default
                        % tensor orientation is with respect to the purple
                        % (1) soliton.
                        newmatcell = matcell_notrot;
                    else
                        newsolitondirection = i;
                        [exx_rot,exy_rot,eyx_rot,eyy_rot,gxy_rot] = obj.rotateStrainTensor(newsolitondirection);
                        newmatcell = {exx_rot,exy_rot,eyx_rot,eyy_rot,gxy_rot,obj.fixed_body_rotation,obj.principal_strain_1,obj.principal_strain_diff};
                        legendcell = {'exx rotated','exy reconstituted','eyx reconstituted','eyy rotated','gxy rotated','fixed body rotation','principal strain max component','maximum (principal) shear'};
                    end
                    
                    for j = 1:nmats
                        this_new_mat = newmatcell{j};
                        these_vals = this_new_mat(SP_only_masks{i});
                        storage_cell{j,i} = these_vals;
                    end
                end
                
                % Having extracted all values from the individual soliton
                % walls, now process into stats. Storage mat extraction
                % functions in the same way as before.
                for j = 1:nmats
                    all_values = vertcat(storage_cell{j,:});
                    storage_mat(j,:) = [mean(all_values),median(all_values),std(all_values),max(all_values),min(all_values)];
                end
                names = legendcell';
                values = storage_mat(:,1);
                stds = storage_mat(:,3);
                npixels = nnz(SP_only_masks{1}) + nnz(SP_only_masks{2}) + nnz(SP_only_masks{3});
                
                % Note: standard error calculation is going to have to be
                % fixed for all directions setting, because original code
                % has recomputation of values on a conncomp basis.
                
            else
                % We may still need to rotate the tensor if we are
                % calculating against an SP other than 1
                %
                % Modification on 07/27/2020: because SP rotation now takes
                % place within the calculatStrain03() method, there is no
                % need to rotate the tensor here. The only thing that is
                % needed is to get the masks along the particular SP
                % direction of interest, and we can either average
                % everything else out or try to back calculate the exact
                % results from the number of components. 
                newmatcell = matcell_notrot;
%                 if SPdirection == 1
%                     % Assumption for the time being is that the default
%                     % tensor orientation is with respect to the purple
%                     % (1) soliton.
%                     newmatcell = matcell_notrot;
%                 else
%                     [exx_rot,exy_rot,eyx_rot,eyy_rot,gxy_rot] = obj.rotateStrainTensor(SPdirection);
%                     newmatcell = {exx_rot,exy_rot,eyx_rot,eyy_rot,gxy_rot,obj.fixed_body_rotation,obj.principal_strain_1,obj.principal_strain_diff};
%                     legendcell = {'exx rotated','exy reconstituted','eyx reconstituted','eyy rotated','gxy rotated','fixed body rotation','principal strain max component','maximum (principal) shear'};
%                 end
                for j = 1:nmats
                    this_mat = newmatcell{j};
                    these_vals = this_mat(SP_only_masks{1});
                    storage_mat(j,:) = [mean(these_vals),median(these_vals),std(these_vals),max(these_vals),min(these_vals)];
                end
                names = legendcell';
                values = storage_mat(:,1);
                stds = storage_mat(:,3);
                npixels = nnz(SP_only_masks{1});
            end
            
            % Standard error calculation should not occur separately for
            % each saddle point direction, but include data from all masks
            % in one pot (if that is the calculation being run).
            if nargout > 4  % need standard errors, so segment mask into connected components
                % This needs to become more sophisticated for running all
                % tensor directions.
                mask_segmentation = bwconncomp(SP_only_masks{1});
                nconn = mask_segmentation.NumObjects;
                component_storage_mat = zeros(nmats,nconn);  % mean only, different conn comps in cols
                for q = 1:nconn
                    this_mask = false(size(SP_only_masks));
                    this_mask(mask_segmentation.PixelIdxList{q}) = true;
                    for s = 1:nmats
                        % matcell has not changed, can run these
                        % calculations for each matrix separately.
                        this_mat = newmatcell{s};
                        these_vals = this_mat(this_mask);
                        component_storage_mat(s,q) = mean(these_vals);
                    end
                end
                % get std for each type of matrix (each row).
                stdevs_comp = std(component_storage_mat,0,2); % second dimension, default normalization
                stderrs = stdevs_comp / sqrt(nconn);
            end
        end
        
        
        
        
        
        % This is a subtractive mask: whatever is leftover after removing
        % AA and SP regions from the image.
        function untrimmed_AB_mask = makeABmask(obj,SPbuffer_distance,AA_center_type,AAradius_value,useStrainFilterMatchedSP)
            [untrimmed_AAmask,~] = getAACircleBufferMasks(obj,AA_center_type,AAradius_value);
            [untrimmed_SPmask] = makeSaddlePointMask(obj,'all',SPbuffer_distance,useStrainFilterMatchedSP);
            mask = true(size(untrimmed_AAmask));
            untrimmed_AB_mask = mask & ~(untrimmed_SPmask | untrimmed_AAmask);
        end
        
        
        
        
        
        % Originally developed for use in the masking functions, 06/23/2020
        % Note that this does not save the strain maps to any instance
        % variables. So more code will need to be developed to enable
        % plotting of strain maps from tensors in different rotational
        % coordinate systems.
        %
        % Original assumption (06/23/2020): saved strain tensor is
        % oriented to the purple soliton wall. So tensor rotation takes
        % place according to differences in the angles.
        %
        % Modified 07/01/2020: input raw angles here th handle the
        % rotational calibration well. 
        % Angular inputs should be in radians.
        function [exx_rot,exy_rot,eyx_rot,eyy_rot,gxy_rot,newbasis] = rotateStrainTensor(obj,old_tensor_angle,new_tensor_angle,oldbasis)
            % Need to get angle between purple and other soliton.
%             [~,~,~,~,~,~,~,basisrotangle_1] = obj.changeDisplacementBasis2(1);
%             [~,~,~,~,~,~,~,basisrotangle_new] = obj.changeDisplacementBasis2(new_soliton_id);

%             if basisrotangle_new < basisrotangle_1
%                 basisrotangle_new = basisrotangle_new + pi;
%             end
            
            % No constraints on CW or CCW direction of rotation
            angle_diff = new_tensor_angle - old_tensor_angle;
            tensor_rotation_angle = angle_diff;
            tensor_rotation_angle_deg = rad2deg(angle_diff); % for a sanity check.
            fprintf('Rotating strain tensor by %.2f degrees to align with soliton wall.\n',tensor_rotation_angle_deg);
            
            % Rotate the basis so that the axes can be plotted correctly
            rotmat = [cos(tensor_rotation_angle), -sin(tensor_rotation_angle); sin(tensor_rotation_angle), cos(tensor_rotation_angle)];
            newbasis = rotmat*oldbasis;
            
            % Use Colin's tensor rotation code to be fast.
            uP = [cos(tensor_rotation_angle) sin(tensor_rotation_angle)];
            exx_rot = uP(1)^2*obj.exx(:) + 2*uP(1)*uP(2)*obj.gxy(:) + uP(2)^2*obj.eyy(:);
            eyy_rot = uP(1)^2*obj.eyy(:) - 2*uP(1)*uP(2)*obj.gxy(:) + uP(2)^2*obj.exx(:);
            gxy_rot = (uP(1)^2-uP(2)^2)*obj.gxy(:) - uP(1)*uP(2)*(obj.exx(:)-obj.eyy(:));
            % FBR stays constant, invariant under coordinate
            % system changes. Can reconstitute the simple shear
            % terms from FBR the same way in rotated coordinate
            % systems as in the original coordinate system.
            %
            % Note that the way FBR is stored varies depending on what
            % stage of the program's execution we are at. Introduce a new
            % variable to check if we are in rad or deg.
            
            % NOTE: NPK found a bad bug here on 07/04/2020. The plus and
            % minus signs were mixed up, so the simple shear values were
            % reconstituted incorrectly. 
            % This has been fixed. But this whole reconstitution method
            % needs to be vetted better before calculating simple shear
            % strain quantities in multiple soliton walls. We need to PROVE
            % that it works.
            if strcmp(obj.FBR_type,'radians')
%                 exy_rot = gxy_rot + obj.fixed_body_rotation(:);
%                 eyx_rot = gxy_rot - obj.fixed_body_rotation(:);
                exy_rot = gxy_rot - obj.fixed_body_rotation(:);
                eyx_rot = gxy_rot + obj.fixed_body_rotation(:);
            elseif strcmp(obj.FBR_type,'degrees')
%                 exy_rot = gxy_rot + deg2rad(obj.fixed_body_rotation(:));
%                 eyx_rot = gxy_rot - deg2rad(obj.fixed_body_rotation(:));
                exy_rot = gxy_rot - deg2rad(obj.fixed_body_rotation(:));
                eyx_rot = gxy_rot + deg2rad(obj.fixed_body_rotation(:));
            end
            
            exx_rot = reshape(exx_rot,size(obj.exx));
            exy_rot = reshape(exy_rot,size(obj.exx));
            eyx_rot = reshape(eyx_rot,size(obj.exx));
            eyy_rot = reshape(eyy_rot,size(obj.exx));
            gxy_rot = reshape(gxy_rot,size(obj.exx));
        end
        
        
        
        
        
        
        function [values,stds,npixels,names,stderrs,nconn] = getStrainInABBufferMask(obj,SPbufferdist,AAremovalradius,AAcentertype,useStrainFilterMatchedSP,strain_crop_range)
            if nargin < 6
                strain_crop_range = [];
            end
            
            untrimmed_AB_mask = obj.makeABmask(SPbufferdist,AAcentertype,AAremovalradius,useStrainFilterMatchedSP);
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            [xspace,yspace] = meshgrid(xbase,ybase);
            [rs,cs] = size(xspace);
            [tr,tc] = size(obj.exx);
            crop_value = (rs-tr)/2;
            
            
            
            
            if ~isempty(strain_crop_range)
                trimmed_AB_mask = trimArray(untrimmed_AB_mask,obj.trim_value);
                trimmed_AB_mask = trimmed_AB_mask(strain_crop_range{1},strain_crop_range{2});
            else
                trimmed_AB_mask = trimArray(untrimmed_AB_mask,crop_value); 
            end
            
            TEST_FIGURES = true;
            if TEST_FIGURES
                figure;
                imagesc(obj.xaxis,obj.yaxis,trimmed_AB_mask);set(gca,'yDir','normal');axis equal;
                title('Trimmed AB mask');
            end
            
            % compute the value for each quantity
            matcell = {obj.exx,obj.exy,obj.eyx,obj.eyy,obj.gxy,obj.fixed_body_rotation,obj.principal_strain_1};
            legendcell = {'exx','exy','eyx','eyy','gxy','fixed body rotation','principal strain max component'};
            nmats = numel(matcell);  % what we have stored and of interest at the moment
            storage_mat = zeros(nmats,5);  % mean, median, std, max, min
            for j = 1:nmats
                this_mat = matcell{j};
                these_vals = this_mat(trimmed_AB_mask);
                storage_mat(j,:) = [mean(these_vals),median(these_vals),std(these_vals),max(these_vals),min(these_vals)];
            end
            names = legendcell';
            values = storage_mat(:,1);
            stds = storage_mat(:,3);
            npixels = nnz(trimmed_AB_mask);
            
            
            if nargout > 4  % need standard errors, so segment mask into connected components
                mask_segmentation = bwconncomp(trimmed_AB_mask);
                nconn = mask_segmentation.NumObjects;
                component_storage_mat = zeros(nmats,nconn);  % mean only, different conn comps in cols
                for q = 1:nconn
                    this_mask = false(size(trimmed_AB_mask));
                    this_mask(mask_segmentation.PixelIdxList{q}) = true;
                    for s = 1:nmats
                        % matcell has not changed, can run these
                        % calculations for each matrix separately.
                        this_mat = matcell{s};
                        these_vals = this_mat(this_mask);
                        component_storage_mat(s,q) = mean(these_vals);
                    end
                end
                % get std for each type of matrix (each row).
                stdevs_comp = std(component_storage_mat,0,2); % second dimension, default normalization
                stderrs = stdevs_comp / sqrt(nconn);
            end
        end
        
        
        
        
        
        function setStrainFilterMatchedSPRegistration(obj,filterstruct,AA_type,AA_buffer_radius,optional_SP_tolerance)
            
            [AA_untrimmed_mask,~] = getAACircleBufferMasks(obj,AA_type,AA_buffer_radius);
            
            SP_hsv = [0,1,1;
                0.33,1,1;
                0.66,1,1];
            AB_hsv = [0,0,1;
                0.33,0,1;
                0.66,0,1];
            % Note that AA_lininds_mask is actually a full mask now.
%             AA_lininds_mask = obj.getAACircleMask(radius_multiplier,optional_AAcircle_min_radius);
            xdisp = obj.annealed_dfield(:,:,1);
            ydisp = obj.annealed_dfield(:,:,2);
            [ xdisp_filt,ydisp_filt ] = filterDisplacement( xdisp,ydisp,filterstruct,false,obj );
            extended_zone_disps = [xdisp_filt(:),ydisp_filt(:)];
            [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( extended_zone_disps );
            xdisp_rz = reshape(reduced_zone_disps(:,1),size(xdisp));
            ydisp_rz = reshape(reduced_zone_disps(:,2),size(ydisp));
            displacement_field = cat(3,xdisp_rz,ydisp_rz);
            [fighn,~] = makeOutsideCustomDisplacementColorPlot(obj,xdisp_filt,ydisp_filt);
            plot_flag = 1;
            upsample_flag = 0;
            SP_id = 1;
            [ SP_line1 ] = registerUniqueSPLines( displacement_field, AA_untrimmed_mask, SP_id, plot_flag, fighn, upsample_flag, optional_SP_tolerance );
            SP_id = 2;
            [ SP_line2 ] = registerUniqueSPLines( displacement_field, AA_untrimmed_mask, SP_id, plot_flag, fighn, upsample_flag, optional_SP_tolerance );
            SP_id = 3;
            [ SP_line3 ] = registerUniqueSPLines( displacement_field, AA_untrimmed_mask, SP_id, plot_flag, fighn, upsample_flag, optional_SP_tolerance );
            obj.SP_registration_unique_strainfiltermatched = {SP_line1,SP_line2,SP_line3}; 
            obj.SP_registration_unique_filter = filterstruct;
            obj.SP_registration_unique_dfield = displacement_field;
        end
        
        
        
        
        
        % Compares the original one used to generate the annealing and the
        % new one produced by the function
        function compareSPRegistrations(obj)
            xdisp = obj.SP_registration_unique_dfield(:,:,1);
            ydisp = obj.SP_registration_unique_dfield(:,:,2);
            [figh,~] = obj.makeOutsideCustomDisplacementColorPlot(xdisp,ydisp);
            hold on 
            SPlines = obj.SP_registration_unique_strainfiltermatched{1} | obj.SP_registration_unique_strainfiltermatched{2} | obj.SP_registration_unique_strainfiltermatched{3};
            h = imagesc(obj.xaxis,obj.yaxis,SPlines);
            h.AlphaData = 0.5;
            title('Post strain filter saddle point registration');
            [figh,~] = obj.makeOutsideCustomDisplacementColorPlot(xdisp,ydisp);
            hold on;
            h2 = imagesc(obj.xaxis,obj.yaxis,full(obj.soliton_walls_merged));
            h2.AlphaData = 0.5;
            title('Pre-annealing saddle point registration');
        end
        
        
        
        
        
        
        
        % This function builds masks on the basis of the triangulation and
        % computes the strain. Intended for DS3 and DS6, which need
        % customized treatment
        function [AAaverages,AAstds,AAmeantwistangles,AAstdtwistangles] = getTriangulationStrainTrends(obj,AAradius_calc)
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            [xspace,yspace] = meshgrid(xbase,ybase);
            [rs,cs] = size(xspace);
            [tr,tc] = size(obj.exx);
            crop_value = (rs-tr)/2;
            
            
            % Build AA masks
            centers = (obj.AAtriangulation.Points-1)*obj.scan_stepsize;  % because all of the AA centers are originally in pixels, stored this way by DT
            % Now anything built off of "centers" will index correctly to
            % the triangulation vertices.
            ncent = size(centers,1);
            % The -1 acccounts that the axes are zero indexed while
            % the pixels are 1 indexed.
            tfstor = false(rs,cs);
            AAmasks = cell(ncent,1);
            for j = 1:ncent
                this_cent = centers(j,:);
                [tf] = isInCircle(xspace,yspace,this_cent(1),this_cent(2),AAradius_calc);
                AAmasks{j} = tf;
            end
            
            % Assign an average twist angle for each AA vertex. Take
            % average of the surrounding six (or fewer triangles).
            AAaveraged_twists = zeros(ncent,1);
            AAstdtwistangles = zeros(ncent,1);
            for i = 1:ncent
                triangle_inds = obj.AAtriangulation.vertexAttachments(i);
                % Get Moire angles corresponding to these inds
                if isempty(triangle_inds{:})
                    AAaveraged_twists(i) = nan;
                    AAstdtwistangles(i) = nan;
                else
                    AAaveraged_twists(i) = mean(obj.triangle_moire_angle(triangle_inds{:}));
                    AAstdtwistangles(i) = std(obj.triangle_moire_angle(triangle_inds{:}));
                end
            end
            
            % Compute quantities in the AA regions
            % Apply the masks to the different components of the strain
            % field.
            matcell = {obj.exx,obj.exy,obj.eyx,obj.eyy,obj.gxy,obj.fixed_body_rotation,obj.principal_strain_1};
            legendcell = {'exx','exy','eyx','eyy','gxy','fixed_body_rotation','principal strain max'};
            nmats = numel(matcell);  % what we have stored and of interest at the moment
            av_storage_mat = zeros(ncent,nmats);
            stdev_storage_mat = zeros(ncent,nmats);
            for i = 1:ncent
                for j = 1:nmats
                    this_mat = matcell{j};
                    this_AAmask = trimArray(AAmasks{i},crop_value);
                    these_vals = this_mat(this_AAmask);
                    av_storage_mat(i,j) = mean(these_vals);
                    stdev_storage_mat(i,j) = std(these_vals);
                end
            end
            AAaverages = av_storage_mat;
            AAstds = stdev_storage_mat;
            AAmeantwistangles = AAaveraged_twists;
        end
        
        
        
        
        
        
        % Function for probing potential connections between uniaxial
        % strain and 2D strain maps
        function mean_linecut_exx = getStrainMapLineCut(obj,target_cell)
            % Plot volumetric strain elements and make selections on top of
            % those.
            load('demo_colormap.mat');
            figure;
            imagesc(trimArray(obj.xaxis,10),trimArray(obj.yaxis,10),obj.exx); axis equal; set(gca,'yDir','normal'); 
            scaleColorMap(cmap,0);
            xlabel('x');
            ylabel('y');
            
%             disp('Please click the start and the end of the line cut.');
%             [x,y] = ginput(2);
            
            disp('Please select the left and right endpoints of the SP region to analyze, parallel along the direction of the soliton wall.');
            coords = ginput(2);
            fc = gcf;
            
            while true
                r = input('Please enter the distance in each direction the fitting region may go: ');
%                 if linedrawn
%                     close(fc);
%                     [fc,axc] = genfig();
%                 end
                figure(fc);
%                 axes(axc);
                Linep = coords(2,:) - coords(1,:);
                parallel_angle = atan(Linep(2)/Linep(1));
                perp_angle = pi/2 + parallel_angle;
                hold on
                inc = [r*cos(perp_angle),r*sin(perp_angle)];
                bcoords1 = coords(1,:) + inc;
                bcoords2 = coords(1,:) - inc;
                bcoords3 = coords(2,:) + inc;
                bcoords4 = coords(2,:) - inc;
                
                l1 = line([bcoords1(1),bcoords2(1)],[bcoords1(2),bcoords2(2)],'Color','g','LineWidth',2.5);
                l2 = line([bcoords2(1),bcoords4(1)],[bcoords2(2),bcoords4(2)],'Color','g','LineWidth',2.5);
                l3 = line([bcoords4(1),bcoords3(1)],[bcoords4(2),bcoords3(2)],'Color','g','LineWidth',2.5);
                l4 = line([bcoords3(1),bcoords1(1)],[bcoords3(2),bcoords1(2)],'Color','g','LineWidth',2.5);
                boundary_coords = vertcat(bcoords1,bcoords2,bcoords3,bcoords4);
                linedrawn = true;
                
                yn = input('Is the displayed ROI defined satisfactorily?');
                outercontinue = 0;
                if ~yn
                    yn = input('Do you want to (0) modify the distance of the rectangular ROI, or (1) redefine the points of the ROI?');
                    if yn
                        outercontinue = 1;
                        break
                    else
                        outercontinue = 0;
                    end
                    
                else
                    % save the roi figure
%                     if saveplot_flag
%                         fh = gcf;
%                         titlename = sprintf('ROI for %dth 1D strain calculation',count);
%                         title(titlename);
%                         savename1 = sprintf('1DstrainROI_%d.fig',count);
%                         savename2 = sprintf('1DstrainROI_%d.png',count);
%                         currentd = pwd;
%                         cd(obj.saved_plots_folderpath);
%                         if ~exist(obj.saved_plots_foldername,'dir')
%                             mkdir(obj.saved_plots_foldername);
%                         end
%                         cd(obj.saved_plots_foldername);
%                         savefig(fh,savename1);
%                         saveas(fh,savename2);
%                         cd(currentd);
%                     end
                    break
                end
            end
            if outercontinue % meaning we want to go back to the definition of the ROI region
%                 continue
            end
            
            
            % Get the lineslice values.
            % Get a linecut average over the soliton
            linecut = @(base_t,N) [linspace(0,r*2,N)'.*cos(perp_angle),linspace(0,r*2,N)'.*sin(perp_angle)] + repmat((1-base_t)*bcoords2 + base_t*bcoords4,N,1);
            
            Ncuts = 50;
            NinCut = 100;
            linecutexx = zeros(NinCut,Ncuts);
            basecoords = linspace(0,1,Ncuts);
            obj.reformAxes();
            xbase = obj.xaxis;
            ybase = obj.yaxis;
            [xspace,yspace] = meshgrid(xbase,ybase);
            
            % for this particular dataset 3
            xspace = trimArray(xspace,10);
            yspace = trimArray(yspace,10);
            
            for i = 1:Ncuts
                this_linecut = linecut(basecoords(i),NinCut);
                %     dispxq = interp2(dfieldx,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
                %     dispyq = interp2(dfieldy,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
                %                     dispxq = interp2(dfieldx_relative,this_linecut(:,1),this_linecut(:,2));
                %                     dispyq = interp2(dfieldy_relative,this_linecut(:,1),this_linecut(:,2));
                exxq = interp2(xspace,yspace,obj.exx,this_linecut(:,1),this_linecut(:,2));
%                 dispyq = interp2(xspace,yspace,dfieldy_relative,this_linecut(:,1),this_linecut(:,2));
                this_line_disps = [exxq];
                linecutexx(:,i) = this_line_disps;
                % save the starting and ending coordinates for later
                % use.
                if i == floor(Ncuts/2)
                    lower_AB_RS_point = this_linecut(1,:);
                    higher_AB_RS_point = this_linecut(end,:);
                end
            end
            mean_linecut_exx = mean(linecutexx,2);
            fh3 = figure;
            plot(1:NinCut,mean_linecut_exx','b-o');
        end
        
        
        
        
        
        
        % Initially, this function will only plot Gaussian ellipse fits,
        % but in the future it may do more than that
        %
        % Stored coordinates are in pixels, not nm.
        function plotAllAAFits(obj,outlier_removal)
            if nargin < 2
                outlier_removal = [];
            end
            xbase = 1:obj.datacube_size(4);
            ybase = 1:obj.datacube_size(3);
            [xspace,yspace] = meshgrid(xbase,ybase);
            space = [xspace(:),yspace(:)];
            ncoords = size(obj.AAgaussianfit_raw_parameters,1);
            net_predvals = zeros(size(space,1),1);
            for i = 1:ncoords
                if nnz(i == outlier_removal) > 0
                    continue
                end
                these_params = obj.AAgaussianfit_raw_parameters(i,1:7); % the eighth entry is the RMSRs.
                x0 = these_params(1);
                y0 = these_params(2);
                sigmax = these_params(3);
                sigmay = these_params(4);
                A = these_params(5); 
                B = these_params(6);
                C = these_params(7);
                [ predvals ] = ellipticGaussianPredfun( space,x0,y0,sigmax,sigmay, A, B, C );
                net_predvals = net_predvals + predvals;
            end
            
            net_predvals = reshape(net_predvals,size(xspace));
            figure;
            imagesc(obj.xaxis,obj.yaxis,-net_predvals); axis equal; set(gca,'yDir','normal');
            colormap(fire); colorbar;
            
        end
        
        
        
        
        
        
        % This function adapated from filter_DSC in the BilayerGraphene.m
        % class.
        % Helper Function for assigning pseudostacking regions
        % Returns a binary mask showing which pixels meet the assigned
        % displacement vector criteria.
        function mask = getPsuedostackingMask(obj,type,r,filterstruct,use_annealed,croprange,manual_disp_input,manual_crop_cell)
            if nargin < 7
                manual_disp_input = [];
            end
            if nargin < 8
                manual_crop_cell = [];
            end
            if use_annealed == 0
                DSC_field = obj.DSC_fit_storage;
                uX = DSC_field(:,:,1);
                uY = DSC_field(:,:,2);
                if croprange > 0
                    uX = trimArray(uX,croprange);
                    uY = trimArray(uY,croprange);
                end
                plotflag = false;
                [ uX,uY ] = filterDisplacement( uX,uY,filterstruct,plotflag,obj );
                
            elseif use_annealed == 1
                DSC_field = obj.annealed_dfield;
                uX = DSC_field(:,:,1);
                uY = DSC_field(:,:,2);
                if croprange > 0
                    uX = trimArray(uX,croprange);
                    uY = trimArray(uY,croprange);
                end
                plotflag = false;
                [ uX,uY ] = filterDisplacement( uX,uY,filterstruct,plotflag,obj );
                
            elseif use_annealed == 2
                DSC_field = obj.dfield_filtered_for_strain*10;
                % The *10 is necessary becuase the values are stored in nm,
                % but we need them in Angstrom.
                uX = DSC_field(:,:,1);
                uY = DSC_field(:,:,2);
                % No cropping or filtering because the strain calculation
                % function has already done this.
            elseif use_annealed == 3
                uX = manual_disp_input(:,:,1);
                uY = manual_disp_input(:,:,2);
            else
                error('Invalid option for "use_annealed"');
            end
               
            % For options 1 and 2, still need to perform the same
            % unwrapping and reconversion to reduced zone.
            if use_annealed > 0
                uXl = uX(:);
                uYl = uY(:);
                [ rzd ] = extendedZoneDisp2ReducedZoneDisp( [uXl,uYl] );
                uX = reshape(rzd(:,1),size(uX));
                uY = reshape(rzd(:,2),size(uY));
            end
            
            
            t = 60/180*pi;
            rotmat = [cos(t) sin(t); -sin(t) cos(t)];
            [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
%             v1 = [0;obj.a/sqrt(3)];
%             v2 = rotmat*v1;
            if nargin < 3
                r = 0.3;
            end
            
            if ~isempty(manual_crop_cell)
                uX = uX(manual_crop_cell{1},manual_crop_cell{2});
                uY = uY(manual_crop_cell{1},manual_crop_cell{2});
            end
            
            switch type
                case 'AA'
                    if strcmp(r,'modified wigner-seitz')
                        r = 0.71;  % For the purposes of AA registration, 
                        % anything with amplitude <= 0.71 is called AA.
                        mask = isInCircle(uX,uY,0,0,r);
                    else
                        mask = isInCircle(uX,uY,0,0,r);
                    end
                case 'AB'
                    if strcmp(r,'modified wigner-seitz')
                        rAA = 0.71;
                        AAmask = isInCircle(uX,uY,0,0,rAA);
                        [ scores ] = SPscore( [uX(:),uY(:)] );
                        reformed_scores = reshape(scores,size(uX));
                        is_AB_angle = reformed_scores <= 0.5;
                        mask = is_AB_angle & ~AAmask;
                    else
                        h(1,:) = v1;
                        h(2,:) = v2;
                        h(3,:) = v2-v1;
                        h(4,:) = -v1;
                        h(5,:) = -v2;
                        h(6,:) = v1-v2;
                        [tf1] = isInCircle(uX,uY,h(1,1),h(1,2),r);
                        [tf2] = isInCircle(uX,uY,h(2,1),h(2,2),r);
                        [tf3] = isInCircle(uX,uY,h(3,1),h(3,2),r);
                        [tf4] = isInCircle(uX,uY,h(4,1),h(4,2),r);
                        [tf5] = isInCircle(uX,uY,h(5,1),h(5,2),r);
                        [tf6] = isInCircle(uX,uY,h(6,1),h(6,2),r);
                        mask = tf1 | tf2 | tf3 | tf4 | tf5 | tf6;
                    end
                case 'SP'
                    if strcmp(r,'modified wigner-seitz')
                        rAA = 0.71;
                        AAmask = isInCircle(uX,uY,0,0,rAA);
                        [ scores ] = SPscore( [uX(:),uY(:)] );
                        reformed_scores = reshape(scores,size(uX));
                        is_SP_angle = reformed_scores > 0.5;
                        mask = is_SP_angle & ~AAmask;
                    else
                        h(1,:) = (v1+v2)/2;
                        h(2,:) = -(v1+v2)/2;
                        h(3,:) = (2*v1-v2)/2;
                        h(4,:) = -(2*v1-v2)/2;
                        h(5,:) = (2*v2-v1)/2;
                        h(6,:) = -(2*v2-v1)/2;
                        [tf1] = isInCircle(uX,uY,h(1,1),h(1,2),r);
                        [tf2] = isInCircle(uX,uY,h(2,1),h(2,2),r);
                        [tf3] = isInCircle(uX,uY,h(3,1),h(3,2),r);
                        [tf4] = isInCircle(uX,uY,h(4,1),h(4,2),r);
                        [tf5] = isInCircle(uX,uY,h(5,1),h(5,2),r);
                        [tf6] = isInCircle(uX,uY,h(6,1),h(6,2),r);
                        mask = tf1 | tf2 | tf3 | tf4 | tf5 | tf6;
                    end
                case 'AA_outer_irevised'
                    router = 0.71;  % For the purposes of AA registration, 
                        % anything with amplitude <= 0.71 is called AA.
                    rinner = 0.35;
                    innermask = isInCircle(uX,uY,0,0,rinner);
                    outermask = isInCircle(uX,uY,0,0,router);
                    mask = outermask & ~innermask;
                case 'AA_inner_irevised'
%                     router = 0.71;  % For the purposes of AA registration, 
                        % anything with amplitude <= 0.71 is called AA.
                    rinner = 0.35;
                    innermask = isInCircle(uX,uY,0,0,rinner);
%                     outermask = isInCircle(uX,uY,0,0,router);
                    mask = innermask;
                case 'AB_irevised'
                    rAA = 0.71;
                    AAmask = isInCircle(uX,uY,0,0,rAA);
                    [ scores ] = SPscore( [uX(:),uY(:)] );
                    reformed_scores = reshape(scores,size(uX));
                    is_AB_angle = reformed_scores <= 0.25;
                    mask = is_AB_angle & ~AAmask;
                case 'SP_irevised'
                    rAA = 0.71;
                    AAmask = isInCircle(uX,uY,0,0,rAA);
                    [ scores ] = SPscore( [uX(:),uY(:)] );
                    reformed_scores = reshape(scores,size(uX));
                    is_SP_angle = reformed_scores > 0.75;
                    mask = is_SP_angle & ~AAmask;
                case 'ABSP_transition_irevised'
                    rAA = 0.71;
                    AAmask = isInCircle(uX,uY,0,0,rAA);
                    [ scores ] = SPscore( [uX(:),uY(:)] );
                    reformed_scores = reshape(scores,size(uX));
                    is_transition_angle = (reformed_scores <= 0.75) & (reformed_scores > 0.25);
                    mask = is_transition_angle & ~AAmask;
            end
        end
        
        
        
        
        
        
        % Note on 07/11/2020: r should have the value 'modified wigner-seitz'
        % in order to use the new partitioning criterion.
        %
        % Modified by NPK on 01/19/2021 in order to enable intermediate
        % stacking ranges to be analyzed in order to address reviewer
        % comments. 
        % In this case, set r = 'intermediate'
        function [AA_perc,AB_perc,SP_perc,pseudostacking_mat] = getPseudostackingAssignments(obj,make_plots,r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell)
            if isempty(r)
                r = 0.3;
            end
            if nargin < 7
                dispfield = [];
            end
            if nargin < 8
                manual_crop_cell = [];
            end
            if strcmp(r,'intermediate')
                AAoutermask = obj.getPsuedostackingMask('AA_outer_irevised',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                ABmask = obj.getPsuedostackingMask('AB_irevised',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                SPmask = obj.getPsuedostackingMask('SP_irevised',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                AAinnermask = obj.getPsuedostackingMask('AA_inner_irevised',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                ABSPtransitionmask = obj.getPsuedostackingMask('ABSP_transition_irevised',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                
                stacking_struct.AAoutermask = AAoutermask;
                stacking_struct.ABmask = ABmask;
                stacking_struct.SPmask = SPmask;
                stacking_struct.AAinnermask = AAinnermask;
                stacking_struct.ABSPtransitionmask = ABSPtransitionmask;
                
                plotmask = zeros(size(AAoutermask));
                plotmask(AAinnermask) = 1;
                plotmask(AAoutermask) = 2;
                plotmask(SPmask) = 3;
                plotmask(ABmask) = 5;
                plotmask(ABSPtransitionmask) = 4;
                
                % Oh, I guess we don't need this.
%                 AA_perc = stacking_struct; % just to get the correct number of arguments out. 
                if make_plots
                    % No explicit cropping or filtering because the strain calculation
                    % function has already done this; use those axes.
                    %                 xax = trimArray(obj.xaxis,croprange);
                    %                 yax = trimArray(obj.yaxis,croprange);
                    if use_annealed > 0
                        xb = obj.xbase_for_strainmaps;
                        yb = obj.ybase_for_strainmaps;
                    else
                        xb = obj.xaxis;
                        yb = obj.xaxis;
                    end
                    figure;
                    imagesc(xb,yb,plotmask); set(gca,'yDir','normal'); axis equal; colorbar; colormap(fire);
                    title('Pseudostacking assignments: 0 = none, 1 = AAinner, 2 = AAouter, 3 = pure SP, 4 = pure AB, 5 = AB/SP transition');
                end
                n_total = numel(AAoutermask);
                
                stacking_struct.AAouterperc = nnz(AAoutermask)/n_total*100;  
                stacking_struct.ABperc = nnz(ABmask)/n_total*100;  
                stacking_struct.SPperc = nnz(SPmask)/n_total*100;  
                stacking_struct.AAinnerperc = nnz(AAinnermask)/n_total*100;
                stacking_struct.ABSPtransitionperc = nnz(ABSPtransitionmask)/n_total*100;  
                
                AA_perc = stacking_struct;
                AB_perc = 'see AA perc struct';
                SP_perc = 'see AA perc struct';
                pseudostacking_mat = plotmask;
                
            else
                AAmask = obj.getPsuedostackingMask('AA',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                ABmask = obj.getPsuedostackingMask('AB',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                SPmask = obj.getPsuedostackingMask('SP',r,filterstruct,use_annealed,croprange,dispfield,manual_crop_cell);
                plotmask = zeros(size(AAmask));
                plotmask(AAmask) = 1;
                plotmask(SPmask) = 2;
                plotmask(ABmask) = 3;
                
                if make_plots
                    % No explicit cropping or filtering because the strain calculation
                    % function has already done this; use those axes.
                    %                 xax = trimArray(obj.xaxis,croprange);
                    %                 yax = trimArray(obj.yaxis,croprange);
                    if use_annealed > 0
                        xb = obj.xbase_for_strainmaps;
                        yb = obj.ybase_for_strainmaps;
                    else
                        xb = obj.xaxis;
                        yb = obj.xaxis;
                    end
                    figure;
                    imagesc(xb,yb,plotmask); set(gca,'yDir','normal'); axis equal; colorbar; colormap(fire);
                    title('Pseudostacking assignments: 0 = none, 1 = AA, 2 = SP, 3 = AB');
                end
                n_total = numel(AAmask);
                AA_perc = nnz(AAmask)/n_total*100;
                AB_perc = nnz(ABmask)/n_total*100;
                SP_perc = nnz(SPmask)/n_total*100;
                pseudostacking_mat = plotmask;
            end
        end
        
        
        
        
        
        % This function is solely for the purpose of determining how the
        % real diffraction pattern is rotated relative to the theoretical
        % diffraction pattern.
        function getDiffractionPatternRotation(obj,method)
            warning('on','all');
            switch method
                case 'disk registration'
                    error('This method not yet supported');
                case 'integration disks'
                    warning('This method may not be optimally accurate.');
                    if isempty(obj.com_coords)  % These are calculated on the basis of the integration disks.
                        obj.getCOMbeamCenter();
                    end
                    % center all of the disks
                    centered_disks = obj.disk_centers - obj.com_coords;
                    % compute the angles
%                     angles = atan2(-centered_disks(:,2),centered_disks(:,1)); 
                    % Note that these angles reflect the flipped y axis on
                    % the DP plots. Roll with this for now -- it's embedded
                    % in the disk coords.
                    
                    % redone so it is the actual angles and actual expected
                    % angles.
                    angles = atan2(centered_disks(:,2),centered_disks(:,1));
                    % angles from horizontal
                    expected_angles = -[pi/2,pi/6,-pi/6,-pi/2,-5*pi/6,5*pi/6,2*pi/3,pi/3,0,-pi/3,-2*pi/3,pi]';
                    % now get the angle differences, accounting for
                    % possible mod issues
                    expected_angles_mod = mod(expected_angles,2*pi);
                    actual_angles_mod1 = mod(angles,2*pi);
                    actual_angles_mod2 = mod(angles,2*pi) + 2*pi;
                    actual_angles_mod3 = mod(angles,2*pi) - 2*pi;
                    diff1 = actual_angles_mod1 - expected_angles_mod;
                    diff2 = actual_angles_mod2 - expected_angles_mod;
                    diff3 = actual_angles_mod3 - expected_angles_mod;
                    % Needs to be actual minus expected because if actual
                    % is larger, then the disk rotation has been positive.
                    combdiffs = [diff1,diff2,diff3];
                    % Want the smallest magnitude angle difference, but
                    % then need the sign again for the angle itself
                    [~,bestidxs] = min(abs(combdiffs),[],2);
                    bestdiffs = zeros(12,1);
                    for i = 1:12
                        bestdiffs(i) = combdiffs(i,bestidxs(i));
                    end
                    rotation_angle = mean(bestdiffs,'omitnan');
                    obj.diffraction_theta_deg = rad2deg(rotation_angle);
                otherwise 
                    error('Please select a valid diffraction pattern rotation estimatino method.');
            end
             
        end
        
        
        
        
        
        % This function takes into account all rotational factors to
        % provide a first-principles rotation of the x-axis so that the
        % derivatives will match the displacement field. SP rotation is
        % returned as a bonus.
        function [realspacerotation_deg,SProtation_deg] = getRealSpaceRotationAngle(obj)
            if isempty(obj.scan_direction_theta_deg)
                error('Please enter m4.scan_direction_theta_deg');
            end
            if isempty(obj.diffraction_theta_deg)
                error('Please run m4.getDiffractionPatternRotation()');
            end
            
            % rotation in degrees
            theta_q = obj.diffraction_theta_deg;
            theta_s = obj.scan_direction_theta_deg;
            realspacerotation_deg = -60 - theta_q - 19 + theta_s;
            % get angle aligned with purple SP
            [~,~,~,~,~,~,~,SProtation_rad] = obj.changeDisplacementBasis2(1);
            SProtation_deg = rad2deg(SProtation_rad);
        end
        
        
        
        
        
        
        
        % Using the same basic approach as 1D strain mapping in terms of
        % the window. Two major exceptions: (1) automatically generated
        % windows based on connected components of SP after AA mask
        % segmentation, and (2) no Sine-Gordon fit.
        %
        % Recommended value for pixel_size_cutoff is around 5. It's the
        % number of SP wall registered pixels that have to be in an SP
        % domain for it to get used.
        function [distcoords,av_scores,distance,endpoint1,endpoint2,rangemean,rangemeanstderr] = getSPwidths(obj,filterstruct,cut_number,cut_spacing_nm,AA_center_type,radius_value,image_boundary_buffer,pixel_size_cutoff,window_dims_nm,cutoff_score)
            TESTFIGS = 1;
            % First attempt is going to be without using the annealed
            % displacement field, to try to keep maximal fidelity to the
            % observed experimental data. This is the approach taken in the
            % AA radius fitting, too.
            dfieldx = obj.annealed_dfield(:,:,1);  % because we are referencing this against the v1
            dfieldy = obj.annealed_dfield(:,:,2); % change this momentarily
            [ dfieldx_extended,dfieldy_extended ] = filterDisplacement( dfieldx,dfieldy,filterstruct,0,obj );
            if nargin < 5
                cut_number = [];
            end
%           If we decide to use the annealed dfield for this calculation, will need these lines of code.             
            [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( [dfieldx(:),dfieldy(:)] );
            dfieldx = reshape(reduced_zone_disps(:,1),size(dfieldx));
            dfieldy = reshape(reduced_zone_disps(:,2),size(dfieldy));
            genfig = @() obj.makeOutsideCustomDisplacementColorPlot(dfieldx,dfieldy);

            % Generate AA masks to chop out of the saddle points. 
            [untrimmed_mask,trimmed_mask] = getAACircleBufferMasks(obj,AA_center_type,radius_value);
            % Show the remaining regions
            if TESTFIGS > 0
                figure; genfig(); hold on
                h = imagesc(obj.xaxis([1,end]),obj.yaxis([1,end]),untrimmed_mask);
                h.AlphaData = 0.5;
            end
            
            % Chop chop
            % For a single run of the function, only one type of SP will be
            % analyzed to avoid this being a mess like some of the other
            % functions here. 
            % Later, add functionality to choose which, but for now do all
            % at once.
            spreg = full(obj.soliton_walls_merged);
            [rsp,csp] = size(spreg);
            [rmask,cmask] = size(untrimmed_mask);
            assert(rsp == rmask);
            assert(csp == cmask);
            spreg(untrimmed_mask) = 0;
            if TESTFIGS > 0
                figure; genfig(); hold on
                h = imagesc(obj.xaxis([1,end]),obj.yaxis([1,end]),spreg);
                h.AlphaData = 0.5;
            end
            
            % Obtain the connected components of the SP network, and filter
            % for distance to the edge of the image.
            spreg_seg = bwconncomp(spreg,8);
            ncomps = spreg_seg.NumObjects;
            sp_to_use_sub = cell(0,1);
            sp_to_use_lin = cell(0,1);
            rlb = image_boundary_buffer;
            clb = image_boundary_buffer;
            rub = rsp - image_boundary_buffer;
            cub = csp - image_boundary_buffer;
            % Likely an off-by-one thing here, but it isn't important at
            % all.
            for i = 1:ncomps
                these_pixels_lin = spreg_seg.PixelIdxList{i};
                [these_pixels_r,these_pixels_c] = ind2sub(size(spreg),these_pixels_lin);
                c1 = any(these_pixels_r < rlb);
                c2 = any(these_pixels_r > rub);
                c3 = any(these_pixels_c < clb);
                c4 = any(these_pixels_c > cub);
                % Clean by removing all SPs with five pixels or fewer (not
                % going to be able to get an accurate slope quantification)
                c5 = numel(these_pixels_lin) <= pixel_size_cutoff;
                if c1 || c2 || c3 || c4 || c5
                    % don't include this to the list of regions to use.
                else
                     sp_to_use_sub{end+1} = [these_pixels_r,these_pixels_c];
                     sp_to_use_lin{end+1} = these_pixels_lin;
                end
            end
            n_clean = length(sp_to_use_sub);
            % Test by visualizing the collected SPs after the edge clean.
            if TESTFIGS > 0
                netsp = zeros(size(spreg));
                for j = 1:n_clean
                     netsp(sp_to_use_lin{j}) = 1;
                end
                figure; genfig(); hold on
                h = imagesc(obj.xaxis([1,end]),obj.yaxis([1,end]),netsp);
                h.AlphaData = 0.5;
            end
            
            % For each connected region, we are now going to ascertain the
            % slope and draw perpendicular averaging lines. SP character
            % data will be stored for each, then we will produce an average
            % and also a cutoff at the end. 
            %
            % Key choice: are we interpolating on HSV character or on xy
            % coordinates here?
            %
            % First try: interpolate on xy and use the annealed dfield.
            % Rationale: this will help normalize differences over fitting
            % region stability and avoid the "avoided crossing" bias, which 
            % would be an artifact rather than a real effect.
            [xspace,yspace] = meshgrid(obj.xaxis,obj.yaxis);
            wwidth = window_dims_nm(1);
            wlength = window_dims_nm(2);
            cutN = wlength/cut_spacing_nm + 1;  % linspace N needs to be fenceposts
            score_storage = zeros(cutN,n_clean);
            distcoords = (-wlength/2:cut_spacing_nm:wlength/2)';
            for i = 1:n_clean
                this_wall = sp_to_use_lin{i}; 
                xcom = mean(xspace(this_wall));
                ycom = mean(yspace(this_wall));
                xcoords_centered = xspace(this_wall) - xcom;
                ycoords_centered = yspace(this_wall) - ycom;
                % These will be linear so we can get the angle
                % directly. Could also do a least squares fit but this
                % might be fine.
                
                % rise over run better for these small SPs
                m = (ycoords_centered(end)-ycoords_centered(1))/(xcoords_centered(end)-xcoords_centered(1));
                % get the unit vector and angle if we need them.
                if (m ~= inf && m ~= -inf)
                    s = [1/sqrt(1+m^2),m/sqrt(1+m^2)];
                else
                    s = [0,1];
                end
                % get the 90 degree normal vector too
                rotmat = [0,-1;1,0];
                s2 = (rotmat*s')';
                angle = atan(s(2)/s(1));  % The extra rotation needs to make this -pi/2
                
                % Build the interpolation window
                
                edgep1 = [xcom,ycom] + wwidth/2*s;
                edgep2 = [xcom,ycom] - wwidth/2*s;
                incvec = wlength/2*s2;
                cor1 = edgep1-incvec;
                cor2 = edgep2-incvec;
                cor3 = edgep2+incvec;
                cor4 = edgep1+incvec;
                if TESTFIGS > 1
                    figure; genfig(); hold on;
                    coords = vertcat(cor1,cor2,cor3,cor4);
                    h = patch(coords(:,1),coords(:,2),'r');
                    h.FaceAlpha = 0.5;
                end
                
                % Within this window, we will now do interpolated slices
%                 % Starting point is the code from Muller's 1D strain
%                 % mapping.
%                 dfieldx_relative = dfieldx;
%                 dfieldy_relative = dfieldy;
%                 % Get a linecut average over the soliton
%                 perpangle = angle;
%                 bcoords2 = cor1;
%                 bcoords4 = cor2;  % I think this is the correct correspondence
%                 
                % produce x and y coords of the linecut
                % believe this is flipped from org
                % Solve for N giving desired stepsize
                
%                 dist = sum((cor1-cor4).^2).^0.5;
%                 dist/
                newlinecut = @(base_t,N) [linspace(cor1(1),cor4(1),cutN)',linspace(cor1(2),cor4(2),cutN)'] + [base_t*(cor2(1)-cor1(1)),base_t*(cor2(2)-cor1(2))];
%                 linecut = @(base_t,N) [linspace(0,r*2,N)'.*cos(perp_angle),linspace(0,r*2,N)'.*sin(perp_angle)] + repmat((1-base_t)*bcoords2 + base_t*bcoords4,N,1);
%                 if isempty(cut_number)
%                     Ncuts = 50;
%                 else
%                     Ncuts = cut_number;
%                 end
%                 NinCut = 100;
                linecutdisps = zeros(cutN,2,cut_number);
%                 basecoords = linspace(0,1,Ncuts);
%                 obj.reformAxes();
%                 xbase = obj.xaxis;
%                 ybase = obj.yaxis;
%                 [xspace,yspace] = meshgrid(xbase,ybase);
                for q = 1:cut_number
%                     this_linecut = linecut(basecoords(i),NinCut);
                    t = q/cut_number;
                    [this_linecut] = newlinecut(t,cut_spacing_nm);
                    %     dispxq = interp2(dfieldx,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
                    %     dispyq = interp2(dfieldy,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
%                     dispxq = interp2(dfieldx_relative,this_linecut(:,1),this_linecut(:,2));
%                     dispyq = interp2(dfieldy_relative,this_linecut(:,1),this_linecut(:,2));
                    dispxq = interp2(xspace,yspace,dfieldx_extended,this_linecut(:,1),this_linecut(:,2));
                    dispyq = interp2(xspace,yspace,dfieldy_extended,this_linecut(:,1),this_linecut(:,2));
                    this_line_disps = [dispxq,dispyq];
                    linecutdisps(:,:,q) = this_line_disps;
%                     % save the starting and ending coordinates for later
%                     % use.
%                     if i == floor(Ncuts/2)
%                         lower_AB_RS_point = this_linecut(1,:);
%                         higher_AB_RS_point = this_linecut(end,:);
%                     end

                    
%%% TODO: for each linecutdisp, save distance coord from the SP COM line,
%%% convert to reduced zone, measure angles, assign SP character, plot,
%%% decide if convolutional alignment is needed
% Looks like the variance here is large enough on a single line cut on a
% ccount of pixel noise that we will need to average after all.
                end
                mean_linecut_disps = mean(linecutdisps,3);
                [ reduced_zone_linedisps ] = extendedZoneDisp2ReducedZoneDisp( [mean_linecut_disps(:,1),mean_linecut_disps(:,2)] );
                [ scores ] = SPscore( reduced_zone_linedisps );    
                score_storage(:,i) = scores;
            end
            
            % For a first method, naively average without alignment. This
            % actually should be pretty good if the SP registration is
            % good. Alignment could run into problems on which SP to
            % convolve first, e.g.
            
            av_scores = mean(score_storage,2);
            figure;
            plot(distcoords,av_scores,'-o');
            xlabel('Distance from SP center (nm)');
            ylabel('SP score (reduced displacement angle)');
            title('Averaged SP line profile');
            
            % For a given SP cutoff (say, 0.5), interpolate a lot and find
            % the one closest to zero.
            interpN = 1000;
            interp_distcoords = linspace(distcoords(1),distcoords(end),interpN);
            much_interpolated = interp1(distcoords,av_scores,interp_distcoords);
            segment1 = much_interpolated(1:round(interpN/2));
            segment2 = much_interpolated((1+round(interpN/2)):end);
            [~,idx1] = min(abs(segment1 - cutoff_score));
            [~,idx2] = min(abs(segment2 - cutoff_score));
            endpoint1 = interp_distcoords(idx1);
            endpoint2 = interp_distcoords(idx2 + round(interpN/2));
            distance = abs(endpoint2 - endpoint1);
            
            % NPK 07/11/2020: for uncertainty quantification, try doing
            % this a different way: calculate the 0.5 cutoff for every
            % single SP used and see what the standard error is. 
            end1stor = zeros(n_clean,1);
            end2stor = zeros(n_clean,1);
            rangestor = zeros(n_clean,1);
            for j = 1:n_clean
                this_profile = score_storage(:,j);
                interpN = 1000;
                interp_distcoords = linspace(distcoords(1),distcoords(end),interpN);
                much_interpolated = interp1(distcoords,this_profile,interp_distcoords);
                segment1 = much_interpolated(1:round(interpN/2));
                segment2 = much_interpolated((1+round(interpN/2)):end);
                [~,idx1] = min(abs(segment1 - cutoff_score));
                [~,idx2] = min(abs(segment2 - cutoff_score));
                this_endpoint1 = interp_distcoords(idx1);
                this_endpoint2 = interp_distcoords(idx2 + round(interpN/2));
                this_distance = abs(this_endpoint2 - this_endpoint1);
                end1stor(j) = this_endpoint1;
                end2stor(j) = this_endpoint2;
                rangestor(j) = this_distance;
            end
            
            % Get average values and stdevs.
            end1mean = mean(end1stor);
            end2mean = mean(end2stor);
            rangemean = mean(rangestor);  % Probably going to be smaller
            rangemeanstderr = std(rangestor)/sqrt(n_clean);
            
            % Sanity check plots
            figure;
            histogram(end1stor);
            xlabel('Endpoint 1 values (nm)'); ylabel('counts'); title('Endpoint 1 location');
            figure;
            histogram(end2stor);
            xlabel('Endpoint 2 values (nm)'); ylabel('counts'); title('Endpoint 2 location');
            figure;
            histogram(rangestor);
            xlabel('Range values (nm)'); ylabel('counts'); title('Mean of ranges');
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Data manipulation methods %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
        
        
        
        
        
        
        function plotSingleDP(obj,idx1,idx2,make_histograms_flag,remove_hBN,remove_Graphene,remove_Beamstop,correct_DP)
            if nargin < 2
                rand1 = randi(obj.datacube_size(3));
                rand2 = randi(obj.datacube_size(4));
            else
                rand1 = idx1;
                rand2 = idx2;
            end
            if nargin < 4
                make_histograms_flag = false;
            end
            if nargin < 5
                remove_hBN = false;
                remove_Graphene = false;
                remove_Beamstop = false;
                correct_DP = false;
            end
            DP = obj.singleLoad(rand1,rand2);
            
            if correct_DP
                subtract_elastic_flag = 1;
                radial_elastic_correction_flag = 0;
                interpolation_correction_flag = 1;
                DP = obj.correctDP(DP,subtract_elastic_flag,radial_elastic_correction_flag,interpolation_correction_flag);
            end
            if remove_hBN
                DP(obj.hBN_mask) = 0;
            end
            if remove_Graphene
                DP(obj.graphene_mask) = 0;
            end
            if remove_Beamstop
                DP(obj.beamstop_mask) = 0;
            end
            
%             figure
            plotDP(DP,1);
            
            figure
            surf(DP); shading flat
            
            if make_histograms_flag
                [all_disk_mask,inner_disk_mask,outer_disk_mask,masks_stack] = obj.getMasksForDiskIntegration();
                for i = 1:12
                    these_counts = DP(masks_stack(:,:,i));
                    edges = -0.25:0.5:max(double(these_counts));
                    figure
                    histogram(these_counts,edges);
                    title(sprintf('Detector counts for %dth disk, DP [%d,%d]',i,rand1,rand2));
                end
            end
        end
        
        
        
        function plotAveragedDP(obj,exponent,loadnum)
            if ~isempty(obj.averaged_DP) && nargin < 3
                dataav = obj.averaged_DP;
                exponentu = exponent;
            else
                if nargin < 2
                    i = 1;
                    exponentu = 1;
                else
                    i = loadnum;
                    exponentu = exponent;
                end
                dataslice = obj.partialLoad(i);
                dataav = squeeze(mean(mean(dataslice,4),3));
                obj.averaged_DP = dataav;
            end
            plotDP(dataav,exponentu);
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
        
        function darkFieldCallback(obj,position,plothandle,plothandle2)
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
            
            % Add the blinking pattern.
            xdisps = obj.DSC_fit_storage(:,:,1);
            ydisps = obj.DSC_fit_storage(:,:,2);
            xavdisp = mean(xdisps(boolmat));
            yavdisp = mean(ydisps(boolmat));
            axes(plothandle2);
            bound_handling_flag = 1;
            scaling_constant = 1;
            [ pred_vals ] = trigFittingFunctions( [xavdisp,yavdisp], scaling_constant, bound_handling_flag );
            
            RADIUS = 0.1;
            [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
            positions = [v1;v2;v2-v1;-v1;-v2;v1-v2;...
                2*v1-v2;v1+v2;2*v2-v1;v2-2*v1;-v1-v2;v1-2*v2];
            n = 500;
            xbase = linspace(-3,3,n);
            ybase = linspace(-3,3,n);
            [xspace,yspace] = meshgrid(xbase,ybase);
            r = 0.2;
            plotmat_randimine = zeros(size(xspace));
            
            for i = 1:12
                [tf] = isInCircle(xspace,yspace,positions(i,1),positions(i,2),r);
                plotmat_randimine(tf) = pred_vals(i);
            end
            
            colormap(fire);
            imagesc(plotmat_randimine);
            shading flat
            colorbar
            caxis([0,1]);
            set(plothandle2,'ydir','normal');
            % xlim([0,500]);
            % ylim([0,500]);
            axis square
        end
        
        
        
        
        function inverseDarkFieldCallback(obj,pos,axh,extended_zone_flag)
            persistent lastred
            if isempty(lastred)
                lastred = [];
            end
            axes(axh);
            % % %             ylb = find(obj.yaxis == ceil(pos(1)));
            % % %             yub = find(obj.yaxis == floor(pos(1)+pos(3)));
            % % %             xlb = find(obj.xaxis == ceil(pos(2)));
            % % %             xub  = find(obj.xaxis == floor(pos(2)+pos(4)));
            ylb = ceil(pos(1));
            yub = floor(pos(1)+pos(3));
            xlb = ceil(pos(2));
            xub = floor(pos(2)+pos(4));
            %             DSC1 = obj.DSC_fit_storage(xlb:xub,ylb:yub,1);
            %             DSC2 = obj.DSC_fit_storage(xlb:xub,ylb:yub,2);
            if extended_zone_flag
                DSC1 = obj.annealed_dfield(xlb:xub,ylb:yub,1);
                DSC2 = obj.annealed_dfield(xlb:xub,ylb:yub,2);
            else
                DSC1 = obj.DSC_fit_storage(xlb:xub,ylb:yub,1);
                DSC2 = obj.DSC_fit_storage(xlb:xub,ylb:yub,2);
            end
            meanamp = mean(mean(sqrt((DSC1.^2 + DSC2.^2))));
            title(sprintf('Mean amplitude in ROI: %.2f',meanamp));
            delete(lastred);
            lastred = scatter(DSC1(:),DSC2(:),'filled','r');
        end
        
        
        
        function [ roi_data,mask ] = roiBoilerplate(obj)
            disks_to_use = obj.disk_averages;
            disks_to_use(isnan(disks_to_use)) = 0;
            overlay = sum(disks_to_use,3);
            figure;
            imagesc(overlay);
            title('Overlay of all disks');
            colormap(gray);
            axis equal
            set(gca,'yDir','normal');
            BW = false(obj.datacube_size(3),obj.datacube_size(4));
            while true
                disp('Please select a region of known AB stacking.');
                thisBW = roipoly();
                BW = BW | thisBW;
                yn = input('Do you wish to select another AB region? 1/0');
                if ~yn
                    break
                end
            end
            roi_data = zeros(nnz(BW),12);
            for i = 1:12
                disk_data = disks_to_use(:,:,i);
                roi_data(:,i) = disk_data(BW);
            end
        end
        
        
        
        
        
        % Control the spread of the displacement population.
        function updateNeighboringDisplacementsOrderControl(obj,first_idx)
            these_idxs = obj.updateNeighboringDisplacements(first_idx);
            while ~isempty(these_idxs)
                next_idxs = [];
                disp(these_idxs)
                for i = 1:numel(these_idxs)
                    out_idxs = obj.updateNeighboringDisplacements(these_idxs(i));
                    next_idxs = [next_idxs,out_idxs];
                end
                these_idxs = next_idxs;
            end
        end
        
        
        
        
        % Make this private in the end
        function selected_idxs = updateNeighboringDisplacements(obj,this_idx)
            nAB = size(obj.AB_centroids,1);
            % compute shortest path to all other AB centroids
            this_pos = obj.emitter_pixel_pos(this_idx,:);
            BW_this = false(obj.datacube_size(3:4));
            BW_this(this_pos(1),this_pos(2)) = true;
            [ v1, v2, ~ ] = getDSCBasisVectors();
            w1 = v1 + v2;
            w2 = 2*v2 - v1;
            %             this_idx = first_idx;  % needs refactored
            selected_idxs = [];
            for i = 1:nAB
                if ~obj.emitter_is_set(i)
                    BW_that = false(obj.datacube_size(3:4));
                    that_pos = obj.emitter_pixel_pos(i,:);
                    BW_that(that_pos(1),that_pos(2)) = true;
                    [ path_pixels ] = bwshortestpath( BW_this, BW_that, 1, 1 );
                    % segment on the basis of each Saddle point
                    segmentation_total = 1;
                    crossed_SP_logical = false(1,3);
                    for j = 1:3
                        SPthick = full(obj.soliton_walls_unique{j});
                        SPthick = SPthick | edge(SPthick);
                        path_segemented = path_pixels;
                        path_segemented(SPthick) = false;
                        resstruct = bwconncomp(path_segemented);
                        % Subtract one to account for the inevitable
                        % self-connectivity even in the absence of
                        % segmentation.
                        if resstruct.NumObjects > 1
                            crossed_SP_logical(j) = true;
                        end
                        segmentation_total = segmentation_total + resstruct.NumObjects - 1;
                    end
                    
                    % 06/01/2020: addition here to ensure we cannot
                    % cross the tear mask in the middle.
                    if ~isempty(obj.tear_mask)
                        tear_segmented = path_pixels;
                        tear_segmented(obj.tear_mask) = false;
                        tearresstruct = bwconncomp(tear_segmented);
                        passes = tearresstruct.NumObjects == 1;
                    else
                        passes = true;
                    end
                    
                    display = false;
                    if display
                        switch segmentation_total
                            case 1
                                warning('Path between AB regions was not segmented by SP walls in setEmitterOrientations()! This could indicate a serious problem.');
                            case 2
                                fprintf('Single segmentation: keeping AB centroid %d.\n',i);
                            otherwise
                                fprintf('Multiple segmentations: discarding AB centroid %d.\n',i);
                        end
                        if passes
                            disp('Did not intersect tear mask');
                        else
                            disp('Intersected tear mask');
                        end
                    end
                    
                    % Block for determining the correct displacement of the
                    % new emitter.
                    if segmentation_total == 2 && passes
                        selected_idxs = [selected_idxs,i];
                        SPidx = find(crossed_SP_logical); % Should only have one true value if we made it here.
                        tv = obj.SP_transition_vectors(SPidx,:);
                        % Logic to ascertain if this vector or its negative
                        % should be applied.
                        d = obj.emitter_displacements(this_idx,:);
                        b1_post = (d + tv + v1)';
                        b2_post = (d + tv - v1)';
                        b1_negt = (d - tv + v1)';
                        b2_negt = (d - tv - v1)';
                        A = [w1',w2'];
                        c1post = A\b1_post;
                        c2post = A\b2_post;
                        c1negt = A\b1_negt;
                        c2negt = A\b2_negt;
                        intTOL = 1e-4;
                        postv_test = [all(abs(c1post - round(c1post)) < intTOL), all(abs(c2post - round(c2post)) < intTOL)];
                        negtv_test = [all(abs(c1negt - round(c1negt)) < intTOL), all(abs(c2negt - round(c2negt)) < intTOL)];
                        assert(nnz([postv_test,negtv_test]) == 1);  % or else the derivation has gone wrong.
                        use_pos = any(postv_test);
                        use_neg = any(negtv_test);
                        if use_pos
                            t = tv;
                        elseif use_neg
                            t = -tv;
                        end
                        new_displacement = d + t;
                        obj.emitter_displacements(i,:) = new_displacement;
                        obj.emitter_is_set(i) = true;
                    end
                end
            end
            % Return the next ones that should be run.
            % Recursive logic for completing the rest of the assignments to
            % instance variables will be computed in level 1 helper
            %             selected_idxs
        end
        
        
    end
    
    
    
    
    methods (Static)
        
        
        % Because of the parallel computing implementation
        function [DSCout,residOut,RMSRout,convergence_storageOut] = fitBlinkingInnards(this_disk_averages,trig_prefactor,nan_handle_flag,...
                weight_vector_local,options)
            
            this_disk_averages(isnan(this_disk_averages)) = 0;  % trust the weighting vector has taken care of it.
            % Crudely, start with the [0,0] guess. May need to multistart
            % this if convergence seems suspect.
            %                     initialDSC_guess = [0,0];
            % lsqnonlin syntax -- computing the vector values of the
            % objective function and not the RMSR.
            %             objfcn = @(DSCvector) (this_disk_averages - blinkingFittingFunction(DSCvector)).*weight_vector;
            %                     trig_prefactor = 0.9;
            if options.SpecifyObjectiveGradient
                third_fourth_disk_flag = 0;  % unless we change this
                objfcn = @(DSCvector) extendedTrigFittingFunctionsResidualsWrapper( DSCvector, trig_prefactor, nan_handle_flag, third_fourth_disk_flag, this_disk_averages, weight_vector_local );
            else
                objfcn = @(DSCvector) (this_disk_averages - trigFittingFunctions(DSCvector,trig_prefactor,nan_handle_flag)).*weight_vector_local;
            end
            lb = [-1.24,0];
            trigfun_flag = 1;
            %             lb = [-1.5,-1.5];
            
            ub = [1.24,1.43];  % really a/sqrt(3) and a/2, I believe.
            %             DSC_fit_result = lsqnonlin(objfcn,initialDSC_guess,lb,ub,options);
            
            [ DSC_fit_result, convergence_storageOut ] = extendedMultistartDiskOptimization( objfcn, lb, ub, options, trigfun_flag, 0 );
            DSCout = permute(DSC_fit_result,[3,1,2]);
            this_residuals = objfcn(DSC_fit_result);
            residOut = permute(this_residuals,[3,1,2]);
            RMSRout = rms(objfcn(DSC_fit_result));
        end
        
    end
end

