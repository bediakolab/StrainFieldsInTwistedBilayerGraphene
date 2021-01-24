classdef TwistedBilayerGrapheneAugmented < BilayerGraphene
    % TwistedBilayerGraphene.m
    %
    % recon_struct should have the following fields:
    % .type = 'Tadmor' or 'NPK', depending on whether the old or new
    % reconstruction functions are desired. 
    % 
    % If Tadmor, 
    % .recon_angle gives the AA rotation
    % .recon_distance gives the AA Gaussian FWHM
    %
    % If NPK,
    % .AA_angle gives the AA rotation
    % .AA_distance gives the AA Gaussian sigma 
    % .AB_angle gives the AB rotation. Can specify this as positive or
    % negative, but the result is the same: counterrotation to the AA
    % domains.
    % .AB_buffer gives the distance from the SP walls that the AB should be
    % buffered.
    % .AB_smooth gives the sigma for the Gaussian convolution kernel used
    % to smooth the edges of the triangle.
    %
    % Nathanael Kazmierczak, Apr 2020, Bediako Lab
    
    properties
        moire_angle_deg
        l1_sublat1_idxs
        l1_sublat2_idxs
        l2_sublat1_idxs
        l2_sublat2_idxs
        l1_sublat1_untrimmed
        l1_sublat2_untrimmed
        l2_sublat1_untrimmed
        l2_sublat2_untrimmed
        noModDSCLat
        layer1_fornomod
        raster_scan_centers
    end
    
    methods
        % At some point in the future, add code here to manage the hBN?
%         function obj = TwistedBilayerGrapheneAugmented(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun)
        
        % See above for old constructor syntax. This is the new version.
        % recon_struct will be parsed to figure out which reconstruction
        % function to use. 
        %
        % heterostrain_struct added on 07/02/2020
        function obj = TwistedBilayerGrapheneAugmented(moire_angle_deg,recon_struct,heterostrain_struct,custom_cell_size)
%             obj.name = sprintf('Moire angle = %.2f, Reconstruction angle = %.2f, Reconstruction distance = %.1f',...
%                 moire_angle_deg,recon_angle,recon_distance);
            if nargin < 3
                heterostrain_struct = [];
            end
            if nargin < 4
                custom_cell_size = [];
            end

            obj.name = sprintf('Moire angle = %.2f',moire_angle_deg);
            obj.moire_angle_deg = moire_angle_deg;
            moire_angle_rad = moire_angle_deg*pi/180;
            
            % Lattice structure of graphene
            aG = 2.461;
            cG = 6.708;   % 2 sheets of graphene = cG
            %             basisG = [ ...
            %                 0   1/3 0;
            %                 0   2/3 0;
            %                 1/2 1/6 0;
            %                 1/2 5/6 0];  % shift to A-A alignment at corners
            obj.a = aG;
            obj.c = cG;
            
            basisG_sublat1 = [ ...
                0   1/3 0;
                1/2 5/6 0];
            
            basisG_sublat2 = [ ...
                0   2/3 0;
                1/2 1/6 0];
            
            % Calculate Moire cell length, real Moire angle
            if moire_angle_rad ~= 0
                numCellsEdgeY = round(1 / (2 * sqrt(3) * tan(moire_angle_rad / 2)) - 0.5);
                thetaMoireActual = 2 * atan(1 / (2 * sqrt(3) * (numCellsEdgeY + 0.5)));
                numCellsEdgeY = abs(numCellsEdgeY);  % chirality testing NPK
            else
                defaultmoire = deg2rad(3);
                numCellsEdgeY = round(1 / (2 * sqrt(3) * tan(defaultmoire / 2)) - 0.5);
                thetaMoireActual = 0;
                numCellsEdgeY = abs(numCellsEdgeY);  % chirality testing NPK
            end
            % thetaMoireActual * 180/pi
            % report the Moire angle
            disp(['Actual Moire angle = ' ...
                sprintf('%.04f',thetaMoireActual*180/pi) ...
                ' degrees'])
            
            % Calculate tiling along X direction - should be a real number
            if moire_angle_rad ~= 0
                numCellsEdgeX = sqrt(3) / (2 * tan(thetaMoireActual / 2)) - 0.5;
            else
                numCellsEdgeX = round(sqrt(3) / (2 * tan(deg2rad(3) / 2)) - 0.5);
            end
            numCellsEdgeX = abs(numCellsEdgeX);  % for chirality testing NPK
            
            % Overall unit cell dimensions in X-Y plane
            xLength = aG*sqrt((numCellsEdgeX + 0.5)^2 + (sqrt(3)/2)^2);
            yLength = aG*sqrt(((numCellsEdgeY + 0.5)*sqrt(3))^2 + (1/2)^2);
            if isempty(custom_cell_size)
                cellDimXY = [xLength yLength];
                % Tile initial graphene lattice
                % NPK modification: explicitly separate the two sublattices
                cellDimG = [1 sqrt(3)]*aG;
                numTile  = [numCellsEdgeX numCellsEdgeY];
            else
                % Added principally for making larger unit cells for
                % heterostrain investigations.
                cellDimXY = custom_cell_size;
                cellDimG = [1 sqrt(3)]*aG;
                cell_x_ratio = custom_cell_size(1)/xLength;
                cell_y_ratio = custom_cell_size(2)/yLength;
                numTile  = [numCellsEdgeX*cell_x_ratio, numCellsEdgeY*cell_y_ratio];
            end
            obj.cellDimXY = cellDimXY;
            
           
            % NPK modification 06/08/2020: add unit cells to make no mod
            % strain field calculations better
%             numTile  = [numCellsEdgeX+5 numCellsEdgeY+5];
%             [bb,aa] = meshgrid(-5:(numTile(2)+0),-5:(numTile(1)+1));
            [bb,aa] = meshgrid(-1:(numTile(2)+0),-1:(numTile(1)+1));
            p = [aa(:) bb(:) zeros(numel(aa(:)),1)];
            % separate out the sublattice
            [bInd_sublat1,pInd_sublat2] = meshgrid(1:size(basisG_sublat1,1),1:size(p,1));
            abc_sublat1 = p(pInd_sublat2(:),:) + basisG_sublat1(bInd_sublat1(:),:);
            [bInd_sublat2,pInd_sublat2] = meshgrid(1:size(basisG_sublat2,1),1:size(p,1));
            abc_sublat2 = p(pInd_sublat2(:),:) + basisG_sublat2(bInd_sublat2(:),:);
            % scale to Cartesian coordinates
            xyInit_sublat1 = [abc_sublat1(:,1)*cellDimG(1) abc_sublat1(:,2)*cellDimG(2)];
            xyInit_sublat2 = [abc_sublat2(:,1)*cellDimG(1) abc_sublat2(:,2)*cellDimG(2)];
            
            % Apply heterostrain
            % Need to generate these coordinates first now
            xyG_layer1_sublat1 = xyInit_sublat1;
            xyG_layer1_sublat2 = xyInit_sublat2;
            xyG_layer2_sublat1 = xyInit_sublat1;
            xyG_layer2_sublat2 = xyInit_sublat2;
            if ~isempty(heterostrain_struct)
                % Convention: subject the first layer to strain and not the
                % second. 
                if strcmp(heterostrain_struct.type,'real space')
                    rot = deg2rad(heterostrain_struct.angle);
                    % This is essentially a coordinate axis rotation to
                    % get in the correct frame for applying the
                    % heterostrain, so we can apply the strain along x
                    % and then rotate back.
                    rotmat = [cos(rot),sin(rot);-sin(rot),cos(rot)];
                    rotmat_inv = [cos(rot),-sin(rot);sin(rot),cos(rot)];  % I distrust explicit inverses
                    xyG_layer1_sublat1_rot = (rotmat*xyG_layer1_sublat1')';
                    xyG_layer1_sublat2_rot = (rotmat*xyG_layer1_sublat2')';
                    strain_amount = heterostrain_struct.amount/100;  % the amount is in percent, so we must divide by 100 first.
                    if ~heterostrain_struct.usePoisson
                        % Later on, include the Poisson ratio here and
                        % perform the strain transformation through a
                        % formal strain tensor (so we can have shear
                        % heterostrain and the like).
                        xyG_layer1_sublat1_rot(:,1) = xyG_layer1_sublat1_rot(:,1)*(1+strain_amount);  % distort x coordinates
                        xyG_layer1_sublat2_rot(:,1) = xyG_layer1_sublat2_rot(:,1)*(1+strain_amount);  % distort x coordinates
                    else
                        delta = 0.16; % Poisson ratio for graphene
                        % WARNING: This is in real space and the absolute
                        % magnitudes are probably all wrong, but easiest
                        % for now.
                        transformation_matrix = [1+strain_amount,0;0,1-strain_amount*delta];
                        xyG_layer1_sublat1_rot = (transformation_matrix*xyG_layer1_sublat1_rot')';
                        xyG_layer1_sublat2_rot = (transformation_matrix*xyG_layer1_sublat2_rot')';
                    end
                    % rotate backwards
                    xyG_layer1_sublat1 = (rotmat_inv*xyG_layer1_sublat1_rot')';
                    xyG_layer1_sublat2 = (rotmat_inv*xyG_layer1_sublat2_rot')';
                end
            end
            
            % Generate the 2 lattices by rotating each one +/- half of the Moire angle
            % NPK: recall that these lattices refer to the different layers
            t = thetaMoireActual / 2;
            m1 = [cos(t) -sin(t); sin(t) cos(t)];
            m2 = [cos(t) sin(t); -sin(t) cos(t)];
            
            % There are two sublattices in each layer.
            % Layer 2 gets rotated the opposite way.
            xyG_layer1_sublat1 = xyG_layer1_sublat1 * m1;
            xyG_layer1_sublat2 = xyG_layer1_sublat2 * m1;
            xyG_layer2_sublat1 = xyG_layer2_sublat1 * m2;
            xyG_layer2_sublat2 = xyG_layer2_sublat2 * m2;
            
% %             % Apply reconstruction
% %             % recon_angle is the full reconstruction in this program
% %             if recon_distance ~= 0 || recon_angle ~= 0
% %                 if use_old_reconfun == 1
% %                     [ xyG_layer1_sublat1 ] = gaussRotateLattice( -recon_angle/2, recon_distance, xyG_layer1_sublat1, cellDimXY );
% %                     [ xyG_layer1_sublat2 ] = gaussRotateLattice( -recon_angle/2, recon_distance, xyG_layer1_sublat2, cellDimXY );
% %                     [ xyG_layer2_sublat1 ] = gaussRotateLattice( recon_angle/2, recon_distance, xyG_layer2_sublat1, cellDimXY );
% %                     [ xyG_layer2_sublat2 ] = gaussRotateLattice( recon_angle/2, recon_distance, xyG_layer2_sublat2, cellDimXY );
% %                 elseif use_old_reconfun == 0
% %                     [ xyG_layer1_sublat1 ] = gaussRotateLattice02( -recon_angle/2, recon_distance, xyG_layer1_sublat1, cellDimXY );
% %                     [ xyG_layer1_sublat2 ] = gaussRotateLattice02( -recon_angle/2, recon_distance, xyG_layer1_sublat2, cellDimXY );
% %                     [ xyG_layer2_sublat1 ] = gaussRotateLattice02( recon_angle/2, recon_distance, xyG_layer2_sublat1, cellDimXY );
% %                     [ xyG_layer2_sublat2 ] = gaussRotateLattice02( recon_angle/2, recon_distance, xyG_layer2_sublat2, cellDimXY );
% %                 elseif use_old_reconfun == -1
% %                     [ xyG_layer1_sublat1 ] = gaussRotateLattice03( -recon_angle/2, recon_distance, xyG_layer1_sublat1, cellDimXY );
% %                     [ xyG_layer1_sublat2 ] = gaussRotateLattice03( -recon_angle/2, recon_distance, xyG_layer1_sublat2, cellDimXY );
% %                     [ xyG_layer2_sublat1 ] = gaussRotateLattice03( recon_angle/2, recon_distance, xyG_layer2_sublat1, cellDimXY );
% %                     [ xyG_layer2_sublat2 ] = gaussRotateLattice03( recon_angle/2, recon_distance, xyG_layer2_sublat2, cellDimXY );                  
% %                 end
% %             end

            % Apply new reconstruction protocols, parsed out of the
            % recon_struct.
            % .type = 'Tadmor' or 'NPK', depending on whether the old or new
    % reconstruction functions are desired. 
    % 
    % If Tadmor, 
    % .recon_angle gives the AA rotation
    % .recon_distance gives the AA Gaussian FWHM
    %
    % If NPK,
    % .AA_angle gives the AA rotation
    % .AA_distance gives the AA Gaussian sigma 
    % .AB_angle gives the AB rotation. Can specify this as positive or
    % negative, but the result is the same: counterrotation to the AA
    % domains.
    % .AB_buffer gives the distance from the SP walls that the AB should be
    % buffered.
    % .AB_smooth gives the sigma for the Gaussian convolution kernel used
    % to smooth the edges of the triangle.
            switch recon_struct.type
                case 'Tadmor'
                    AA_angle = recon_struct.recon_angle;
                    AA_distance = recon_struct.recon_distance;
                    [ xyG_layer1_sublat1 ] = gaussRotateLattice03( -AA_angle/2, AA_distance, xyG_layer1_sublat1, cellDimXY );
                    [ xyG_layer1_sublat2 ] = gaussRotateLattice03( -AA_angle/2, AA_distance, xyG_layer1_sublat2, cellDimXY );
                    [ xyG_layer2_sublat1 ] = gaussRotateLattice03( AA_angle/2, AA_distance, xyG_layer2_sublat1, cellDimXY );
                    [ xyG_layer2_sublat2 ] = gaussRotateLattice03( AA_angle/2, AA_distance, xyG_layer2_sublat2, cellDimXY );   
                case 'NPK'
                    AA_angle = recon_struct.AA_angle;
                    AA_distance = recon_struct.AA_distance;
                    AB_angle = recon_struct.AB_angle;
                    AB_buffer = recon_struct.AB_buffer;
                    AB_smooth = recon_struct.AB_smooth;
                    plotrotfield = true;
                    if AB_angle > 0
                        disp('Converting AB_angle to negative value for counterrotation.');
                        AB_angle = -AB_angle;
                    end
                    [ xyG_layer1_sublat1, ~ ] = reconstructionFunctionNPK( -AA_angle/2, AA_distance, -AB_angle/2, AB_buffer, AB_smooth, xyG_layer1_sublat1, cellDimXY, plotrotfield );
                    plotrotfield = false;
                    [ xyG_layer1_sublat2] = reconstructionFunctionNPK( -AA_angle/2, AA_distance, -AB_angle/2, AB_buffer, AB_smooth, xyG_layer1_sublat2, cellDimXY, plotrotfield );
                    plotrotfield = true;
                    [ xyG_layer2_sublat1] = reconstructionFunctionNPK( AA_angle/2, AA_distance, AB_angle/2, AB_buffer, AB_smooth, xyG_layer2_sublat1, cellDimXY, plotrotfield );
                    plotrotfield = false;
                    [ xyG_layer2_sublat2] = reconstructionFunctionNPK( AA_angle/2, AA_distance, AB_angle/2, AB_buffer, AB_smooth, xyG_layer2_sublat2, cellDimXY, plotrotfield );
                case 'NPK_v2'
                    if recon_struct.plotrotfield
                        prf1 = true;
                        prf2 = false;
                    else
                        prf1 = false;
                        prf2 = false;
                    end
                        
                    AA_angle = recon_struct.AA_angle;
                    AA_distance = recon_struct.AA_distance;
                    AB_angle = recon_struct.AB_angle;
                    AB_buffer = recon_struct.AB_buffer;
                    AB_smooth = recon_struct.AB_smooth;
                    gamma = recon_struct.boundary_rectification.gamma;
                    corner_angle_deg = recon_struct.boundary_rectification.corner_angle_deg;
                    plotrotfield = prf1;
                    if AB_angle > 0
                        disp('Converting AB_angle to negative value for counterrotation.');
                        AB_angle = -AB_angle;
                    end
                    [ xyG_layer1_sublat1] = reconstructionFunctionNPK2( -AA_angle/2, AA_distance, -AB_angle/2, AB_buffer, AB_smooth, xyG_layer1_sublat1, cellDimXY, plotrotfield, gamma, corner_angle_deg );
                    plotrotfield = prf2;
                    [ xyG_layer1_sublat2] = reconstructionFunctionNPK2( -AA_angle/2, AA_distance, -AB_angle/2, AB_buffer, AB_smooth, xyG_layer1_sublat2, cellDimXY, plotrotfield, gamma, corner_angle_deg );
                    plotrotfield = prf1;
                    [ xyG_layer2_sublat1] = reconstructionFunctionNPK2( AA_angle/2, AA_distance, AB_angle/2, AB_buffer, AB_smooth, xyG_layer2_sublat1, cellDimXY, plotrotfield, gamma, corner_angle_deg );
                    plotrotfield = prf2;
                    [ xyG_layer2_sublat2] = reconstructionFunctionNPK2( AA_angle/2, AA_distance, AB_angle/2, AB_buffer, AB_smooth, xyG_layer2_sublat2, cellDimXY, plotrotfield, gamma, corner_angle_deg );
                case 'none'
                otherwise
                    error('Please enter a valid recon_struct.type');
            end
            
            
            
            % index the atoms in case we want to later compute the no mod
            % displacement
            xyG_layer1_sublat1_indices = (1:size(xyG_layer1_sublat1,1))';
            xyG_layer1_sublat2_indices = (1:size(xyG_layer1_sublat2,1))';
            xyG_layer2_sublat1_indices = (1:size(xyG_layer2_sublat1,1))';
            xyG_layer2_sublat2_indices = (1:size(xyG_layer2_sublat2,1))';
            
            % Save untrimmed arrays
            obj.l1_sublat1_untrimmed = xyG_layer1_sublat1;
            obj.l1_sublat2_untrimmed = xyG_layer1_sublat2;
            obj.l2_sublat1_untrimmed = xyG_layer2_sublat1;
            obj.l2_sublat2_untrimmed = xyG_layer2_sublat2;
            
            % Remove atoms outside of bounds
            [xyG_layer1_sublat1,xyG_layer1_sublat1_indices] = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer1_sublat1,xyG_layer1_sublat1_indices,cellDimXY);
            [xyG_layer1_sublat2,xyG_layer1_sublat2_indices] = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer1_sublat2,xyG_layer1_sublat2_indices,cellDimXY);
            [xyG_layer2_sublat1,xyG_layer2_sublat1_indices] = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer2_sublat1,xyG_layer2_sublat1_indices,cellDimXY);
            [xyG_layer2_sublat2,xyG_layer2_sublat2_indices] = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer2_sublat2,xyG_layer2_sublat2_indices,cellDimXY);
            
            obj.l1_sublat1 = xyG_layer1_sublat1;
            obj.l1_sublat2 = xyG_layer1_sublat2;
            obj.l2_sublat1 = xyG_layer2_sublat1;
            obj.l2_sublat2 = xyG_layer2_sublat2;
            
            obj.l1_sublat1_idxs = xyG_layer1_sublat1_indices;
            obj.l1_sublat2_idxs = xyG_layer1_sublat2_indices;
            obj.l2_sublat1_idxs = xyG_layer2_sublat1_indices;
            obj.l2_sublat2_idxs = xyG_layer2_sublat2_indices;
            %             del = xyG1(:,1) < 0 ...
            %                 | xyG1(:,1) >= cellDimXY(1) ...
            %                 | xyG1(:,2) < 0 ...
            %                 | xyG1(:,2) >= cellDimXY(2);
            %             xyG1(del,:) = [];
            %             del = xyG2(:,1) < 0 ...
            %                 | xyG2(:,1) >= cellDimXY(1) ...
            %                 | xyG2(:,2) < 0 ...
            %                 | xyG2(:,2) >= cellDimXY(2);
            %             xyG2(del,:) = [];
            
        end
        
        % Method for running the multislice Prism simulations.
        % Will set up a point calculation because this is a trblg, no need
        % to raster.
        function [stack4D,atoms] = simulate(obj,include_hBN_flag)
            [emdSTEM,atoms,cellDim] = simulate@BilayerGraphene(obj,include_hBN_flag);
            disp('Beginning PRISM multislice')
            clear control_struct;
            control_struct.type = 'point';
            emdSTEM = PRISMmultisliceNPK_ExternalControl(emdSTEM,control_struct);
            [stack4D] = PRISMmultiStack(emdSTEM);
            obj.DP = stack4D;  % Just a single diffraction pattern whenever we are doing point calculations.
        end
        
        
        function layer1idxs = getLayer1idxs(obj)
            layer1idxs = vertcat(obj.l1_sublat1_idxs,obj.l1_sublat2_idxs);
        end
        
        function layer2idxs = getLayer2idxs(obj)
            layer2idxs = vertcat(obj.l2_sublat1_idxs,obj.l2_sublat2_idxs);
        end
        
        function layer1_untrimmed = getLayer1_untrimmed(obj)
            layer1_untrimmed = vertcat(obj.l1_sublat1_untrimmed,obj.l1_sublat2_untrimmed);
        end
        
        function layer2_untrimmed = getLayer2_untrimmed(obj)
            layer2_untrimmed = vertcat(obj.l2_sublat1_untrimmed,obj.l2_sublat2_untrimmed);
        end
        
        
        
        % For use in strain calculations. It has to be the specific atom!
        function computeNoModDSCField(obj)
            layer1idxs = obj.getLayer1idxs();
            layer2idxs = obj.getLayer2idxs();
            layer1 = obj.getLayer1_untrimmed();
            layer2 = obj.getLayer2_untrimmed();
%             del = layer1idxs == -1 | layer2idxs == -1;
            layer1c = layer1;
%             layer1c(del,:) = [];
            layer2c = layer2;
%             layer2c(del,:) = [];
            
            obj.noModDSCLat = layer1c - layer2c;
            obj.layer1_fornomod = layer1c;
        end
        
        
        
        function figh = plotNoModDSCField(obj,figh)
            DSC_lat = obj.noModDSCLat;
            tblg_lat1 = obj.layer1_fornomod;
            %             quiver(obj.l1_sublat1(:,1),obj.l1_sublat1(:,2),obj.DSC_sublat1(:,1),obj.DSC_sublat1(:,2));
            %             quiver(obj.l1_sublat2(:,1),obj.l1_sublat2(:,2),obj.DSC_sublat2(:,1),obj.DSC_sublat2(:,2));
            
            DSCamp = sum(DSC_lat.^2,2).^0.5;
            DSCangle = atan2(DSC_lat(:,2),DSC_lat(:,1));
            % figure; contourf(DSCamp,50,'LineStyle','None');
            
            myInterpAmp = scatteredInterpolant(tblg_lat1(:,1),tblg_lat1(:,2),DSCamp,'nearest','nearest');
            myInterpAngle = scatteredInterpolant(tblg_lat1(:,1),tblg_lat1(:,2),DSCangle,'nearest','nearest');
            % make the interpolation points
            npoints = 500;
            xbase = linspace(0,obj.cellDimXY(1),npoints);
            ybase = linspace(0,obj.cellDimXY(2),npoints+1);
            [xspace,yspace] = meshgrid(xbase,ybase);
            interpAmp = myInterpAmp(xspace,yspace);
            interpAngle = myInterpAngle(xspace,yspace);
            
            % Having fixed the DSC boundary case errors by mod padding,
            % this is no longer needed and will make the plots look better.
            %             % Trim three pixels along each edge.
            %             trim_width = 5;
            %             [ interpAmp ] = trimArray( interpAmp,trim_width );
            %             [ interpAngle ] = trimArray( interpAngle,trim_width );
            %             [ xbase ] = trimArray( xbase,trim_width );
            %             [ ybase ] = trimArray( ybase,trim_width );
            
            % Transpose to make this equivalent to Colin's plots.
            figh = figure; contourf(ybase,xbase,interpAmp',50,'LineStyle','None');
            colormap(fire)
            axis equal
            shading interp
            c = colorbar;
            c.Label.String = 'DSC vector magnitude (Angstroms)';
            c.FontSize = 11;
            title(sprintf(strcat(['DSC field amplitude:\n',obj.name])));
            set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('y (Angstroms)');
            ylabel('x (Angstroms)');
            set(gca,'ydir','reverse')
            
            figure; contourf(ybase,xbase,interpAngle',50,'LineStyle','None');
            colormap(hsv)
            axis equal
            shading interp
            c = colorbar;
            c.Label.String = 'DSC vector angle (radians)';
            c.FontSize = 11;
            title(sprintf(strcat(['DSC field angle:\n',obj.name])));
            set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('y (Angstroms)');
            ylabel('x (Angstroms)');
            set(gca,'ydir','reverse')
        end
        
        
        
        
        % Compute strain field according to gradient of the displacement
        % field, as laid out in lecture notes NPK found on 02/01/2020.
        % Please clean this up into proper OOP style at some point.
        function [exx,eyy,exy,eyx] = computeStrainField(obj,remove_mean)
            DSC_field_values = obj.noModDSCLat;
            layer1 = obj.layer1_fornomod;
            % Looks like we need a rectangular grid before we can call
            % gradient().
            % Interpolate both the x and y displacement fields.
            uXInterpolant = scatteredInterpolant(layer1(:,1),layer1(:,2),DSC_field_values(:,1));
            uYInterpolant = scatteredInterpolant(layer1(:,1),layer1(:,2),DSC_field_values(:,2));
            spacing = 0.5;
            xbase = 0:spacing:obj.cellDimXY(1);
            ybase = 0:spacing:obj.cellDimXY(2);
            [obj.strain_xspace,obj.strain_yspace] = meshgrid(xbase,ybase);
            uX = uXInterpolant(obj.strain_xspace,obj.strain_yspace);
            uY = uYInterpolant(obj.strain_xspace,obj.strain_yspace);
            % Now we can compute the gradient with respect to both elements
            % of the displacement field, u(x,y):
            [exx,exy] = gradient(uX,spacing);
            [eyx,eyy] = gradient(uY,spacing);
            
            if remove_mean
                exx = exx - mean(mean(exx));
                eyy = eyy - mean(mean(eyy));
                exy = exy - mean(mean(exy));
                eyx = eyx - mean(mean(eyx));
            end
            filter = false;
            if filter
            filter_cutoff = 0.03;
            exx(exx > filter_cutoff | exx < -filter_cutoff) = 0;
            eyy(eyy > filter_cutoff | eyy < -filter_cutoff) = 0;
            exy(exy > filter_cutoff | exy < -filter_cutoff) = 0;
            eyx(eyx > filter_cutoff | eyx < -filter_cutoff) = 0;
            end
            if remove_mean
                exx = exx - mean(mean(exx));
                eyy = eyy - mean(mean(eyy));
                exy = exy - mean(mean(exy));
                eyx = eyx - mean(mean(eyx));
            end
            obj.exx = exx;
            obj.eyy = eyy;
            obj.exy = exy;
            obj.eyx = eyx;
        end
        
        
        
        function figh = plotStrainField(obj,figh,type)
            if nargin < 2
                figh = figure;
                type = 'all';
            else
                if ~isempty(figh)
                    figure(figh);
                else
                    figh = figure;
                end
            end
            % exx
            %             xbase = obj.strain_xspace(1,:);
            %             ybase = obj.strain_yspace(:,1);
            if strcmp(type,'all') || strcmp(type,'exx')
                % NPK 06/08/2020: changed these from contourf to imagesc
                % because Matlab kept crashing when trying to save the
                % contour plots.
%                 contourf(obj.strain_yspace,obj.strain_xspace,obj.exx*100,50,'LineStyle','None');
               
                yends = [obj.strain_xspace(1,1),obj.strain_xspace(1,end)];
                xends = [obj.strain_yspace(1,1),obj.strain_yspace(end,1)];

                imagesc(xends,yends,obj.exx'*100);
                colormap(fire)
                axis equal
                c = colorbar;
                c.Label.String = 'exx strain %';
                c.FontSize = 11;
                title(sprintf(strcat(['exx strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 700, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'eyy')
%                 contourf(obj.strain_yspace,obj.strain_xspace,obj.eyy*100,50,'LineStyle','None');
                
                imagesc(xends,yends,obj.eyy'*100);
                colormap(fire)
                axis equal
                c = colorbar;
                c.Label.String = 'eyy strain %';
                c.FontSize = 11;
                title(sprintf(strcat(['eyy strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 700, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'exy')
                % Plot the exy strain here (unmodified parameterization)
%                 contourf(obj.strain_yspace,obj.strain_xspace,obj.exy*100,50,'LineStyle','None');
                
                imagesc(xends,yends,obj.exy'*100);
                colormap(fire)
                axis equal
                c = colorbar;
                c.Label.String = 'exy strain %';
                c.FontSize = 11;
                title(sprintf(strcat(['exy strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 700, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
                
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'eyx')
                % Plot the eyx strain here (unmodified parameterization)
%                 contourf(obj.strain_yspace,obj.strain_xspace,obj.eyx*100,50,'LineStyle','None');
                
                imagesc(xends,yends,obj.eyx'*100);
                colormap(fire)
                axis equal
                c = colorbar;
                c.Label.String = 'eyx strain %';
                c.FontSize = 11;
                title(sprintf(strcat(['eyx strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 700, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'gxy')
                % Plot the total shear strain here.
                
                %06/08/2020: NPK added division by half here.
                
                gxy = 0.5*(obj.exy + obj.eyx);
%                 contourf(obj.strain_yspace,obj.strain_xspace,gxy*100,50,'LineStyle','None');
                imagesc(xends,yends,gxy'*100);
                colormap(fire)
                axis equal
                c = colorbar;
                c.Label.String = 'gxy strain %';
                c.FontSize = 11;
                title(sprintf(strcat(['Total shear strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 700, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'theta')
                % Negative sign is because unreconstructed Moire should
                % have positive rotation.
                
                theta = -0.5*(obj.eyx - obj.exy);
%                 contourf(obj.strain_yspace,obj.strain_xspace,rad2deg(theta),50,'LineStyle','None');
                
                imagesc(xends,yends,rad2deg(theta'));
                axis equal
                colormap(fire)
                
                c = colorbar;
                c.Label.String = 'theta (degrees)';
                c.FontSize = 11;
                title(sprintf(strcat(['Fixed-body interlayer rotation:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 700, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            
            
            
        end
        
        
        % NPK added 03/01/2020
        function makeInteractiveDisplacementCircles(obj)
            if isempty(obj.DSC_sublat1) || isempty(obj.noModDSCLat)
                newalgflag = 0;
                obj.computeDSCField(newalgflag);
                obj.computeNoModDSCField();
            end
            figh = figure;
            subplot(1,2,1);
            figh = obj.plotDSCField(figh);
%             obj.DSC
%             damp = (dfield(:,:,1).^2 + dfield(:,:,2).^2).^(1/2);
%             dangle = atan2(dfield(:,:,2),dfield(:,:,1));
%             figure;
%             imagesc(damp); colormap(fire); axis equal; colorbar;
            figs = get(groot,'children');
            close(figs(figs ~= figh));
            figure(figh);
            axh = gca;
            s2 = subplot(1,2,2);
            e = imellipse(axh);
            
            fcn = @(pos) obj.visualizationCallback(pos,s2);
            fcn(e.getPosition());
            id = addNewPositionCallback(e,fcn);
            
        end
        
        
        function visualizationCallback(obj,pos,axh)
            persistent lastscat
            if isempty(lastscat)
                lastscat = [];
            end
            plotFullDisplacementHexagons( axh );
            xmin = pos(1);
            ymin = pos(2);
            xwidth = pos(3);
            ywidth = pos(4);
            x0 = xmin + 1/2*xwidth;
            y0 = ymin + 1/2*ywidth;
            radius = mean([xwidth,ywidth]);
            DSC = obj.noModDSCLat;
            DSC_positions = obj.layer1_fornomod;
            [tf] = isInCircle(DSC_positions(:,1),DSC_positions(:,2),y0,x0,radius);
            DSC_filtered = DSC(tf,:);
            if ~isempty(lastscat)
                delete(lastscat);
            end
            lastscat = scatter(DSC_filtered(:,1),DSC_filtered(:,2),'filled','b');
            
        end
        
        
        
        
        
        % NPK added 04/21/2020 for machine learning collaboration
        % input angle in degrees
        % raster_start_pos in angstrom
        % scan_step in nm
        function [scan_centers,scan_centers_mod] = getRasterScan(obj,raster_start_pos,scan_step,scan_angle,scan_dimensions)
            scan_step = scan_step*10;
            % Build raster scan coordinates
            [xs,ys] = meshgrid(0:(scan_dimensions(1)-1),0:(scan_dimensions(2)-1));
            xs_nm = xs*scan_step;
            ys_nm = ys*scan_step;
            torot = [xs_nm(:),ys_nm(:)];
            ua = deg2rad(scan_angle);
            rotmat = [cos(ua),-sin(ua);sin(ua),cos(ua)];
            rotcoords = torot*rotmat';
%             xs_rotated = reshape(rotcoords(:,1),scan_dimensions);
%             ys_rotated = reshape(rotcoords(:,2),scan_dimensions);
            scan_centers = rotcoords + raster_start_pos;
            obj.raster_scan_centers = scan_centers;
            scan_centers_mod(:,1) = mod(scan_centers(:,1),obj.cellDimXY(1));
            scan_centers_mod(:,2) = mod(scan_centers(:,2),obj.cellDimXY(2));
        end
        
        
        
        
        % input beam_FWHM in nm, converts to A
        function interferometry_dataset = computeBlinkingWeightedAverage(obj,beam_FWHM,max_intensity)
            beam_FWHM = beam_FWHM*10;
            sigma = beam_FWHM/(2*sqrt(2*log(2)));
            cutoff_r = sigma*sqrt(2*log(10)); % half width at 1/10th max
            % Pad the displacement field on each side with replicating
            % displacement vectors 
            displacement_field = obj.getDSClat();
%             org_displacement_field = obj.getDSClat();
            coordinates = obj.getLayer1();
            
            
            pt = coordinates(:,1) > obj.cellDimXY(1) - cutoff_r;
            pb = coordinates(:,1) < cutoff_r;
            pl = coordinates(:,2) > obj.cellDimXY(2) - cutoff_r;
            pr = coordinates(:,2) < cutoff_r;
            % Pad top
            coordinates_padded = coordinates;
            coordinates_padded = vertcat(coordinates_padded, coordinates(pt,:)-[obj.cellDimXY(1),0]);
            displacement_field = vertcat(displacement_field, displacement_field(pt,:));
            % Pad bottom
            coordinates_padded = vertcat(coordinates_padded, coordinates(pb,:)+[obj.cellDimXY(1),0]);
            displacement_field = vertcat(displacement_field, displacement_field(pb,:));
            % Pad left
            coordinates_padded = vertcat(coordinates_padded, coordinates(pl,:)-[0,obj.cellDimXY(2)]);
            displacement_field = vertcat(displacement_field, displacement_field(pl,:));
            % Pad right
            coordinates_padded = vertcat(coordinates_padded, coordinates(pr,:)+[0,obj.cellDimXY(2)]);
            displacement_field = vertcat(displacement_field, displacement_field(pr,:));
            % Pad top left corner
            coordinates_padded = vertcat(coordinates_padded, coordinates(pt&pl,:)-[obj.cellDimXY]);
            displacement_field = vertcat(displacement_field, displacement_field(pt&pl,:));
            % Pad top right corner
            coordinates_padded = vertcat(coordinates_padded, coordinates(pt&pr,:)+[-obj.cellDimXY(1),obj.cellDimXY(2)]);
            displacement_field = vertcat(displacement_field, displacement_field(pt&pr,:));
            % Pad bottom right corner
            coordinates_padded = vertcat(coordinates_padded, coordinates(pb&pr,:)+[obj.cellDimXY]);
            displacement_field = vertcat(displacement_field, displacement_field(pb&pr,:));
            % Pad bottom left corner
            coordinates_padded = vertcat(coordinates_padded, coordinates(pb&pl,:)+[obj.cellDimXY(1),-obj.cellDimXY(2)]);
            displacement_field = vertcat(displacement_field, displacement_field(pb&pl,:));
            
            % Compute beam weighted averages where the beam center is
            % modulo the unit cell dimensions. 
            n_beam_pos = size(obj.raster_scan_centers,1);
            interferometry_dataset = zeros(n_beam_pos,12);
            for i = 1:n_beam_pos
                if mod(i,1000) == 0
                    fprintf('Averaging beam position %d of %d.\n',i,n_beam_pos);
                end
                this_beam_pos = obj.raster_scan_centers(i,:);
                this_beam_pos(1) = mod(this_beam_pos(1),obj.cellDimXY(1));
                this_beam_pos(2) = mod(this_beam_pos(2),obj.cellDimXY(2));
                r = sum((coordinates_padded - this_beam_pos).^2,2).^0.5;
                inds = r < cutoff_r;
                % For displacement vectors located within the cutoff
                % threshold, compute the interferometry vectors.
                to_compute = displacement_field(inds,:);
                radii = r(inds);
                weights = exp(-radii.^2./(2*sigma.^2));
                normweights = weights./sum(sum(weights));
                [r,c] = size(to_compute);
                pred_vals = HamishTrigFittingFunction( to_compute', ones(12,1)*max_intensity );
                pred_vals_rs = reshape(pred_vals,[12,r])';
                averaged_pattern = sum(pred_vals_rs.*normweights,1);
                interferometry_dataset(i,:) = averaged_pattern;
            end
        end
        
        
    end
    
    methods (Static)
        % Refactored for use in the sublattice graphene construction.
        function [atoms,indices] = removeAtomsOutOfCell(atoms,indices,cellDimXY)
            del = atoms(:,1) < 0 ...
                | atoms(:,1) >= cellDimXY(1) ...
                | atoms(:,2) < 0 ...
                | atoms(:,2) >= cellDimXY(2);
            atoms(del,:) = [];
            indices(del,:) = -1;
        end
    end
end

