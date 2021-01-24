classdef TwistedBilayerGraphene < BilayerGraphene
    % TwistedBilayerGraphene.m
    %
    % Modified 07/02/2020 to try to take into account the effects of
    % heterostrain, to gain intuition about this.
    %
    % Nathanael Kazmierczak, Bediako Lab, 2019-2020
    
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
    end
    
    methods
        % At some point in the future, add code here to manage the hBN?
        function obj = TwistedBilayerGraphene(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun)
            if nargin < 4
                use_old_reconfun = 0;
            end
            obj.name = sprintf('Moire angle = %.2f, Reconstruction angle = %.2f, Reconstruction distance = %.1f',...
                moire_angle_deg,recon_angle,recon_distance);
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
            numCellsEdgeY = round(1 / (2 * sqrt(3) * tan(moire_angle_rad / 2)) - 0.5);
            thetaMoireActual = 2 * atan(1 / (2 * sqrt(3) * (numCellsEdgeY + 0.5)));
            % thetaMoireActual * 180/pi
            % report the Moire angle
            disp(['Actual Moire angle = ' ...
                sprintf('%.04f',thetaMoireActual*180/pi) ...
                ' degrees'])
            
            % Calculate tiling along X direction - should be a real number
            numCellsEdgeX = sqrt(3) / (2 * tan(thetaMoireActual / 2)) - 0.5;
            
            % Overall unit cell dimensions in X-Y plane
            xLength = aG*sqrt((numCellsEdgeX + 0.5)^2 + (sqrt(3)/2)^2);
            yLength = aG*sqrt(((numCellsEdgeY + 0.5)*sqrt(3))^2 + (1/2)^2);
            cellDimXY = [xLength yLength];
            obj.cellDimXY = cellDimXY;
            
            % Tile initial graphene lattice
            % NPK modification: explicitly separate the two sublattices
            cellDimG = [1 sqrt(3)]*aG;
            numTile  = [numCellsEdgeX numCellsEdgeY];
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
            
            % Generate the 2 lattices by rotating each one +/- half of the Moire angle
            % NPK: recall that these lattices refer to the different layers
            t = thetaMoireActual / 2;
            m1 = [cos(t) -sin(t); sin(t) cos(t)];
            m2 = [cos(t) sin(t); -sin(t) cos(t)];
            
            % There are two sublattices in each layer.
            % Layer 2 gets rotated the opposite way.
            xyG_layer1_sublat1 = xyInit_sublat1 * m1;
            xyG_layer1_sublat2 = xyInit_sublat2 * m1;
            xyG_layer2_sublat1 = xyInit_sublat1 * m2;
            xyG_layer2_sublat2 = xyInit_sublat2 * m2;
            
            % Apply reconstruction
            % recon_angle is the full reconstruction in this program
            if recon_distance ~= 0 || recon_angle ~= 0
                if use_old_reconfun
                    [ xyG_layer1_sublat1 ] = gaussRotateLattice( -recon_angle/2, recon_distance, xyG_layer1_sublat1, cellDimXY );
                    [ xyG_layer1_sublat2 ] = gaussRotateLattice( -recon_angle/2, recon_distance, xyG_layer1_sublat2, cellDimXY );
                    [ xyG_layer2_sublat1 ] = gaussRotateLattice( recon_angle/2, recon_distance, xyG_layer2_sublat1, cellDimXY );
                    [ xyG_layer2_sublat2 ] = gaussRotateLattice( recon_angle/2, recon_distance, xyG_layer2_sublat2, cellDimXY );
                else
                    % NPK 05/26/2020 update with new reconstruction
                    % function
%                     [ xyG_layer1_sublat1 ] = gaussRotateLattice02( -recon_angle/2, recon_distance, xyG_layer1_sublat1, cellDimXY );
%                     [ xyG_layer1_sublat2 ] = gaussRotateLattice02( -recon_angle/2, recon_distance, xyG_layer1_sublat2, cellDimXY );
%                     [ xyG_layer2_sublat1 ] = gaussRotateLattice02( recon_angle/2, recon_distance, xyG_layer2_sublat1, cellDimXY );
%                     [ xyG_layer2_sublat2 ] = gaussRotateLattice02( recon_angle/2, recon_distance, xyG_layer2_sublat2, cellDimXY );
                    [ xyG_layer1_sublat1 ] = gaussRotateLattice03( -recon_angle/2, recon_distance, xyG_layer1_sublat1, cellDimXY );
                    [ xyG_layer1_sublat2 ] = gaussRotateLattice03( -recon_angle/2, recon_distance, xyG_layer1_sublat2, cellDimXY );
                    [ xyG_layer2_sublat1 ] = gaussRotateLattice03( recon_angle/2, recon_distance, xyG_layer2_sublat1, cellDimXY );
                    [ xyG_layer2_sublat2 ] = gaussRotateLattice03( recon_angle/2, recon_distance, xyG_layer2_sublat2, cellDimXY );
                end
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
            del = layer1idxs == -1 | layer2idxs == -1;
            layer1c = layer1;
            layer1c(del,:) = [];
            layer2c = layer2;
            layer2c(del,:) = [];
            
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
            
            % compute x and y components here too
            xDSC = interpAmp.*cos(interpAngle);
            yDSC = interpAmp.*sin(interpAngle);
            figure; contourf(ybase,xbase,xDSC',50,'LineStyle','None');
            colormap(hsv)
            axis equal
            shading interp
            c = colorbar;
            c.Label.String = 'Angstrom';
            c.FontSize = 11;
            title(sprintf(strcat(['DSC x component:\n',obj.name])));
            set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('y (Angstroms)');
            ylabel('x (Angstroms)');
            set(gca,'ydir','reverse')
            
            figure; contourf(ybase,xbase,yDSC',50,'LineStyle','None');
            colormap(hsv)
            axis equal
            shading interp
            c = colorbar;
            c.Label.String = 'Angstrom';
            c.FontSize = 11;
            title(sprintf(strcat(['DSC y component:\n',obj.name])));
            set(gcf, 'Position',  [0, 0, 500, 800])
            xlabel('y (Angstroms)');
            ylabel('x (Angstroms)');
            set(gca,'ydir','reverse')
        end
        
        
        
        
        % Compute strain field according to gradient of the displacement
        % field, as laid out in lecture notes NPK found on 02/01/2020.
        % Please clean this up into proper OOP style at some point.
        function [exx,eyy,exy,eyx] = computeStrainField(obj,remove_mean)
            DSC_field_values = obj.noModDSCLat();
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
            filter_cutoff = 0.03;
            exx(exx > filter_cutoff | exx < -filter_cutoff) = 0;
            eyy(eyy > filter_cutoff | eyy < -filter_cutoff) = 0;
            exy(exy > filter_cutoff | exy < -filter_cutoff) = 0;
            eyx(eyx > filter_cutoff | eyx < -filter_cutoff) = 0;
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
                contourf(obj.strain_yspace,obj.strain_xspace,obj.exx,50,'LineStyle','None');
                colormap(fire)
                axis equal
                shading interp
                c = colorbar;
                c.Label.String = 'exx';
                c.FontSize = 11;
                title(sprintf(strcat(['exx strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 500, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'eyy')
                contourf(obj.strain_yspace,obj.strain_xspace,obj.eyy,50,'LineStyle','None');
                colormap(fire)
                axis equal
                shading interp
                c = colorbar;
                c.Label.String = 'eyy';
                c.FontSize = 11;
                title(sprintf(strcat(['eyy strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 500, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'exy')
                % Plot the exy strain here (unmodified parameterization)
                contourf(obj.strain_yspace,obj.strain_xspace,obj.exy,50,'LineStyle','None');
                colormap(fire)
                axis equal
                shading interp
                c = colorbar;
                c.Label.String = 'exy';
                c.FontSize = 11;
                title(sprintf(strcat(['exy strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 500, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
                
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'eyx')
                % Plot the eyx strain here (unmodified parameterization)
                contourf(obj.strain_yspace,obj.strain_xspace,obj.eyx,50,'LineStyle','None');
                colormap(fire)
                axis equal
                shading interp
                c = colorbar;
                c.Label.String = 'eyx';
                c.FontSize = 11;
                title(sprintf(strcat(['eyx strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 500, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'gxy')
                % Plot the total shear strain here.
                gxy = obj.exy + obj.eyx;
                contourf(obj.strain_yspace,obj.strain_xspace,gxy,50,'LineStyle','None');
                colormap(fire)
                axis equal
                shading interp
                c = colorbar;
                c.Label.String = 'gxy';
                c.FontSize = 11;
                title(sprintf(strcat(['Total shear strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 500, 800])
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
                set(gca,'ydir','reverse');
            end
            
            if strcmp(type,'all')
                figure;
            end
            if strcmp(type,'all') || strcmp(type,'theta')
                theta = 0.5*(obj.eyx - obj.exy);
                contourf(obj.strain_yspace,obj.strain_xspace,theta,50,'LineStyle','None');
                colormap(fire)
                axis equal
                shading interp
                c = colorbar;
                c.Label.String = 'theta';
                c.FontSize = 11;
                title(sprintf(strcat(['Torsional strain amplitude:\n',obj.name])));
                set(gcf, 'Position',  [0, 0, 500, 800])
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

