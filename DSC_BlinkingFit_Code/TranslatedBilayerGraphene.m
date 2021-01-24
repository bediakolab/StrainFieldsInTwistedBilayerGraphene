classdef TranslatedBilayerGraphene < BilayerGraphene
    % Nathanael Kazmierczak, Dec 20, 2019
    %
    % Class for building unrotated graphene bilayers with a translated
    % offset, enabling electron diffraction simulations on psuedostacked
    % regions. 
    
    properties
        DSC_vector
        DP  % Diffraction pattern from multislice simulation
        averaging_circles  % for getting average intensity of the diffraction disks.
        % should be struct with radii and disk centers
    end
    
    methods
        function obj = TranslatedBilayerGraphene(DSC_vector_or_type,moire_cell_to_match)
            % First, select numeric DSC vector
            obj.a = 2.461;
            obj.c = 6.708;
            if ischar(DSC_vector_or_type) || isstring(DSC_vector_or_type)
                switch DSC_vector_or_type
                    case 'AA'
                        DSC_vector = [0,0];
                    case 'AB'
                        vec1 = [0,obj.a/sqrt(3)];
                        DSC_vector = vec1;
                    case 'SP'
                        t = 60/180*pi;
                        rotmat = [cos(t) sin(t); -sin(t) cos(t)];
                        vec1 = [0,obj.a/sqrt(3)];
                        vec2 = vec1*rotmat';
                        DSC_vector = (vec1+vec2)/2;
                    otherwise
                        error('That string input to the TranslatedBilayerGraphene constructor is not recognized.'); 
                end
                obj.name = sprintf('DSC stacking type: %s',...
                    DSC_vector_or_type);
            else  % presume numeric input
                DSC_vector = DSC_vector_or_type;
                obj.name = sprintf('DSC vector: [%.3f, %.3f]',...
                    DSC_vector_or_type(1),DSC_vector_or_type(2));
            end
            obj.DSC_vector = DSC_vector;
                
            % Generate a tiled hexagonal lattice, as in twisted bilayer
            % construction.
            % The moire cell dimension calculations are still performed to
            % avoid computational artifacts, for now, just in case there is
            % some scaling problem.
            moire_angle_rad = moire_cell_to_match*pi/180;
            
             % Lattice structure of graphene
            aG = 2.461;
            cG = 6.708;   % 2 sheets of graphene = cG
            
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
            % NPK mod: for translation simulations, want it to be an
            % integer multiple of the graphene lattice constant.
            xLength_real = round(xLength/(2*obj.a));
            yLength_real = round(yLength/(obj.a*sqrt(3)));
%             xLength_real = 5;
%             yLength_real = 5;
%             cellDimXY = [xLength yLength];
            cellDimXY = [xLength_real*(2*obj.a) yLength_real*(obj.a*sqrt(3))];
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
            
%             disp('hi');
%             plotTBLG(xyInit_sublat1,xyInit_sublat2,cellDimXY)

            % Separate out the different layers
            xyG_layer1_sublat1 = xyInit_sublat1;
            xyG_layer1_sublat2 = xyInit_sublat2;
            % This needs to be adding the DSC vector because it is defined
            % as moving from the top layer atom to the bottom layer atom. 
            xyG_layer2_sublat1 = xyInit_sublat1 + DSC_vector;
            xyG_layer2_sublat2 = xyInit_sublat2 + DSC_vector;
            
            % Remove atoms outside of bounds
            xyG_layer1_sublat1 = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer1_sublat1,[],cellDimXY);
            xyG_layer1_sublat2 = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer1_sublat2,[],cellDimXY);
            xyG_layer2_sublat1 = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer2_sublat1,[],cellDimXY);
            xyG_layer2_sublat2 = TwistedBilayerGraphene.removeAtomsOutOfCell(xyG_layer2_sublat2,[],cellDimXY);
           
            obj.l1_sublat1 = xyG_layer1_sublat1;
            obj.l1_sublat2 = xyG_layer1_sublat2;
            obj.l2_sublat1 = xyG_layer2_sublat1;
            obj.l2_sublat2 = xyG_layer2_sublat2;
            %  plotTBLG([xyInit_sublat1;xyInit_sublat2],[xyG_layer2_sublat1;xyG_layer2_sublat2],cellDimXY)
            
        end % End constructor
        
        
        % Method for running the multislice Prism simulations.
        % Will set up a point calculation because this is a trblg, no need
        % to raster.
        function [stack4D,atoms] = simulate(obj,include_hBN_flag,extra_spacing_factor)
            [emdSTEM,atoms,cellDim] = simulate@BilayerGraphene(obj,include_hBN_flag,extra_spacing_factor);
            disp('Beginning PRISM multislice')
            clear control_struct;
            control_struct.type = 'point';
            emdSTEM = PRISMmultisliceNPK_ExternalControl(emdSTEM,control_struct);
            [stack4D] = PRISMmultiStack(emdSTEM);
            obj.DP = stack4D;  % Just a single diffraction pattern whenever we are doing point calculations.
        end
        
        
        function setAveragingCircles(obj,disk_centers,radius)
            obj.averaging_circles.disk_centers = disk_centers;
            obj.averaging_circles.radius = radius;
        end
        
        
        function plotAveragingCircles(obj)
            obj.plotDP();
            viscircles(obj.averaging_circles.disk_centers,...
                obj.averaging_circles.radius*ones(size(obj.averaging_circles.disk_centers,1),1));
            xlim([1.725233236151603e+02 4.227390670553936e+02]);
            ylim([1.725233236151603e+02 4.227390670553936e+02]);
            axis equal
        end
        
        
        function av_intensities = getAverageDiskIntensities(obj)
            if isempty(obj.averaging_circles)
                error('Please populate the averaging_circles property using setAveragingCircles prior to invoking getAverageDiskIntensities.');
            end
            
            xbase = 1:size(obj.DP,1);
            ybase = 1:size(obj.DP,2);
            [xspace,yspace] = meshgrid(xbase,ybase);
            
            av_intensities = zeros(size(obj.averaging_circles.disk_centers,1),1);
            for i = 1:size(obj.averaging_circles.disk_centers,1)
                x0 = obj.averaging_circles.disk_centers(i,1);
                y0 = obj.averaging_circles.disk_centers(i,2);
                tf = isInCircle(xspace,yspace,x0,y0,obj.averaging_circles.radius);
                av_intensities(i) = mean(obj.DP(tf));
            end
        end
        
        
        % Plot the lone DP with external code.
        function plotDP(obj,restrict_flag)
            % NPK's hard-coded masking way of getting rid of the central
            % beam
%             addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit');
%             load('12292019_block_simulation_center_beam_mask.mat');
            DP = obj.DP;
%             DP(r) = 0; % Block off the center beam with the mask.
            exponent = 1;
            plotDP(DP,exponent);
            if restrict_flag
                
                xlim([200,400]);
                ylim([180,420]);
                caxis([0,5.5E-6]);
                axis equal
            end
            title(obj.name);
        end
        
    end
    
end

