classdef BilayerGraphene < handle
    % Superclass for TwistedBilayerGraphene and TranslatedBilayerGraphene
    %
    % Nathanael Kazmierczak, Dec 2019
    
    properties
        % The fundamental atomic coordinates. These should be vertical
        % arrays.
        l1_sublat1
        l1_sublat2
        l2_sublat1
        l2_sublat2
        cellDimXY
        name
        DSC_sublat1
        DSC_sublat2
        a % lattice constants
        c 
        exx
        eyy
        exy
        eyx  % strain field values
        strain_xspace
        strain_yspace
    end
    
    methods
        function layer1 = getLayer1(obj)
            layer1 = vertcat(obj.l1_sublat1,obj.l1_sublat2);
        end
        
        function layer2 = getLayer2(obj)
            layer2 = vertcat(obj.l2_sublat1,obj.l2_sublat2);
        end
        
        function DSC_lat = getDSClat(obj)
            if isempty(obj.DSC_sublat1) || isempty(obj.DSC_sublat2)
                warning('At least one DSC sublattice has not been computed. Request to return DSCsublat terminating unsucessfully...');
            else
                DSC_lat = vertcat(obj.DSC_sublat1,obj.DSC_sublat2);
            end
        end
        
        % Function to format and output cell appropriate for the
        % multislice simulations. All of this uses code taken from
        % buildCell02 by Colin Ophus.
        function [atoms,cellDim] = makeBNformatForSimulation(obj)
            
            xyG1 = obj.getLayer1();
            xyG2 = obj.getLayer2();
            
            aBN = 2.504;
            cBN = 6.661;
            basisBN = [ ...
                0   1/3 1/4 5;
                0   2/3 1/4 7;
                1/2 1/6 1/4 7;
                1/2 5/6 1/4 5;
                0   1/3 3/4 7;
                0   2/3 3/4 5;
                1/2 1/6 3/4 5;
                1/2 5/6 3/4 7];
            thetaBN = 14 * pi / 180;
            numBNlayers = 7*1;

            % BN lattice
            uBN = [aBN 0];
            vBN = [0 aBN*sqrt(3)];
            uv0 = [uBN; vBN];
            % Rotate UC vectors
            m = [cos(thetaBN) -sin(thetaBN);
                sin(thetaBN) cos(thetaBN)];
            uvCell = [ ...
                obj.cellDimXY(1) 0;
                0 obj.cellDimXY(2)] * m;
            abTile = round(uv0' \ uvCell')';
            % uProj = abTile(1,1)*cellDimXY(1) + abTile(2,1)*cellDimXY(
            aRange = [min(abTile(:,1)) max(abTile(:,1))] + [-1 1];
            bRange = [min(abTile(:,2)) max(abTile(:,2))] + [-1 1];
            abRangeMax = round(max([abs(aRange) abs(bRange)])*sqrt(3));
            [bb,aa] = meshgrid(-abRangeMax:abRangeMax);
            p = [aa(:) bb(:) zeros(numel(aa),2)];
            [bInd,pInd] = meshgrid(1:size(basisBN,1),1:size(p,1));
            atomsBN = basisBN(bInd(:),:) + p(pInd(:),:);
            atomsBN(:,1:3) = atomsBN(:,1:3).*[aBN aBN*sqrt(3) cBN];
            % Reproject BN into new coordinate system (rational approximate)
            uProj = abTile(1,1)*[aBN 0] + abTile(1,2)*[0 aBN*sqrt(3)];
            vProj = abTile(2,1)*[aBN 0] + abTile(2,2)*[0 aBN*sqrt(3)];
            u0 = [obj.cellDimXY(1) 0];
            v0 = [0 obj.cellDimXY(2)];
            mProj = [uProj; vProj] \ [u0; v0];
            atomsBN(:,1:2) = atomsBN(:,1:2) * mProj;
            % Remove BN atoms outside of cell
            del = atomsBN(:,1) < 0 ...
                | atomsBN(:,1) >= obj.cellDimXY(1) ...
                | atomsBN(:,2) < 0 ...
                | atomsBN(:,2) >= obj.cellDimXY(2);
            atomsBN(del,:) = [];
            
            % Assemble output cell
            cellDimZ = numBNlayers * cBN + 2*obj.c;
            cellDim = [obj.cellDimXY cellDimZ];
            numBN = size(atomsBN,1);
            numG = [size(xyG1,1) size(xyG2,1)];
            atoms = zeros(numBN * numBNlayers + numG(1) + numG(2),4);
            for a0 = 1:numBNlayers
                inds = (1:numBN) + (a0-1)*numBN;
                dz = (a0-1)*cBN;
                atoms(inds,:) = atomsBN ...
                    + [0 0 dz 0];
            end
            
            inds = numBN*numBNlayers + (1:numG(1));
            z = cBN*numBNlayers + obj.c*0.25;
            atoms(inds,1:2) = xyG1;
            atoms(inds,3) = z;
            atoms(inds,4) = 6;
            
            inds = numBN*numBNlayers + numG(1) + (1:numG(2));
            z = cBN*numBNlayers + obj.c*0.75;
            atoms(inds,1:2) = xyG2;
            atoms(inds,3) = z;
            atoms(inds,4) = 6;
        end
        
        
        function [atoms,cellDim] = makeGrapheneForSimulation(obj,extra_spacing_factor)
            xyG1 = obj.getLayer1();
            xyG2 = obj.getLayer2();
            
            % Assemble output cell
            cellDimZ = 2*obj.c;
            cellDim = [obj.cellDimXY cellDimZ];
            numG = [size(xyG1,1) size(xyG2,1)];
            atoms = zeros(numG(1) + numG(2),4);
            
            % Modified NPK on 04/11/2020 to be able to introduce an extra
            % spacing factor in between the atoms.
            inds = 1:numG(1);
            z1 = obj.c*0.25*(1/(1+extra_spacing_factor/2));
            atoms(inds,1:2) = xyG1;
            atoms(inds,3) = z1;
            atoms(inds,4) = 6;
            
            inds = (1:numG(2)) + numG(1);
            z2 = obj.c*0.75*(1+extra_spacing_factor/2);
            atoms(inds,1:2) = xyG2;
            atoms(inds,3) = z2;
            atoms(inds,4) = 6;
        end
        
        
        % Function for running a simulation on a bilayer cell.
        % This is a virtual method -- should not be invoked from
        % superclass. Subclasses will overwrite to make this a
        % point-calculation (trblg) or a raster calculation (tblg) as is
        % appropriate.
        function [emdSTEM,atoms,cellDim] = simulate(obj,include_hBN_flag,extra_spacing_factor)
            addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/Code_from_Colin/BlinkingSimulations/bilayergraphenesimulations');
            if include_hBN_flag
                [atoms,cellDim] = obj.makeBNformatForSimulation();
            else
                [atoms,cellDim] = obj.makeGrapheneForSimulation(extra_spacing_factor);
            end
            emdSTEM = PRISM01(atoms,cellDim);
        end
        
        
        function [emdSTEM,atoms,cellDim] = simulateWithoutHBN(obj,extra_spacing_factor)
            addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/Code_from_Colin/BlinkingSimulations/bilayergraphenesimulations');
            [atoms,cellDim] = obj.makeGrapheneForSimulation(extra_spacing_factor);
            emdSTEM = PRISM01(atoms,cellDim);
        end
        
        
        % Function for computing the DSC lattice
        % See Ishikawa et al. 2016, "Interfacial Atomic Structure of
        % Twisted Few-Layer Graphene", Scientific Reports
        % To avoid outliers, wrap the lattice around within a padding
        % distance.
        %
        % new_alg_flag leads to faster but less accurate computation for
        % small-angle Moire patterns.
        function [DSC_sublat1,DSC_sublat2] = computeDSCField(obj,new_alg_flag)
%             x0 = tblg.l1_sublat1(1,1);
%             y0 = tblg.l1_sublat1(1,2);
            WRAP_MARGIN = 5;
            numlat1 = size(obj.l1_sublat1,1);
            DSC_sublat1 = zeros(size(obj.l1_sublat1));
            l2sl1add1 = obj.l2_sublat1( obj.l2_sublat1(:,1) > (obj.cellDimXY(1) - WRAP_MARGIN),: );
            l2sl1add1 = l2sl1add1 - repmat([obj.cellDimXY(1),0],size(l2sl1add1,1),1);
            l2sl1add2 = obj.l2_sublat1( obj.l2_sublat1(:,1) < WRAP_MARGIN,: );
            l2sl1add2 = l2sl1add2 + repmat([obj.cellDimXY(1),0],size(l2sl1add2,1),1);
            l2_sublat1_ww1 = [obj.l2_sublat1;l2sl1add1;l2sl1add2];
            % having wrapped on X, now this will successfully fill in the
            % corners as well.
            l2sl1add3 = l2_sublat1_ww1( l2_sublat1_ww1(:,2) > (obj.cellDimXY(2) - WRAP_MARGIN),: );
            l2sl1add3 = l2sl1add3 - repmat([0,obj.cellDimXY(2)],size(l2sl1add3,1),1);
            l2sl1add4 = l2_sublat1_ww1( l2_sublat1_ww1(:,2) < WRAP_MARGIN,: );
            l2sl1add4 = l2sl1add4 + repmat([0,obj.cellDimXY(2)],size(l2sl1add4,1),1);
            
            l2_sublat1_withwrap = [l2_sublat1_ww1;l2sl1add3;l2sl1add4];
            
%             figure
%             hold on
%             scatter(obj.l2_sublat1(:,1),obj.l2_sublat1(:,2))
%             scatter(l2sl1add1(:,1),l2sl1add1(:,2))
%             scatter(l2sl1add2(:,1),l2sl1add2(:,2))
%             scatter(l2sl1add3(:,1),l2sl1add3(:,2))
%             scatter(l2sl1add4(:,1),l2sl1add4(:,2))
%             figure
%             scatter(l2_sublat1_withwrap(:,1),l2_sublat1_withwrap(:,2));


            
            
            % Perform the DSC construction process again with the second
            % sublattice.
            l2sl2add1 = obj.l2_sublat2( obj.l2_sublat2(:,1) > (obj.cellDimXY(1) - WRAP_MARGIN),: );
            l2sl2add1 = l2sl2add1 - repmat([obj.cellDimXY(1),0],size(l2sl2add1,1),1);
            l2sl2add2 = obj.l2_sublat2( obj.l2_sublat2(:,1) < WRAP_MARGIN,: );
            l2sl2add2 = l2sl2add2 + repmat([obj.cellDimXY(1),0],size(l2sl2add2,1),1);
            l2_sublat2_ww1 = [obj.l2_sublat2;l2sl2add1;l2sl2add2];
            
            l2sl2add3 = l2_sublat2_ww1( l2_sublat2_ww1(:,2) > (obj.cellDimXY(2) - WRAP_MARGIN),: );
            l2sl2add3 = l2sl2add3 - repmat([0,obj.cellDimXY(2)],size(l2sl2add3,1),1);
            l2sl2add4 = l2_sublat2_ww1( l2_sublat2_ww1(:,2) < WRAP_MARGIN,: );
            l2sl2add4 = l2sl2add4 + repmat([0,obj.cellDimXY(2)],size(l2sl2add4,1),1);
            
            l2_sublat2_withwrap = [l2_sublat2_ww1;l2sl2add3;l2sl2add4];
%             
%             figure
%             hold on
%             scatter(obj.l2_sublat2(:,1),obj.l2_sublat2(:,2))
%             scatter(l2sl2add1(:,1),l2sl2add1(:,2))
%             scatter(l2sl2add2(:,1),l2sl2add2(:,2))
%             scatter(l2sl2add3(:,1),l2sl2add3(:,2))
%             scatter(l2sl2add4(:,1),l2sl2add4(:,2))
%             figure
%             scatter(l2_sublat2_withwrap(:,1),l2_sublat2_withwrap(:,2));
            
            
            DIV_INCREMENT = 10;
%             if obj.cellDimXY(1) < DIV_INCREMENT
%                 ALG = 'old';
%             else
%                 ALG = 'new';
%             end
            if new_alg_flag
                ALG = 'new';
            else
                ALG = 'old';
            end
            PADDING2 = WRAP_MARGIN + 0.1;
            switch ALG
                case 'new'
                    
                    incs = (obj.cellDimXY+2*PADDING2)./DIV_INCREMENT;
                    xcrit = (obj.cellDimXY(1)+PADDING2):-incs(1):-PADDING2;
                    ycrit = (obj.cellDimXY(2)+PADDING2):-incs(2):-PADDING2;
                    
                    for i = 1:4  % 4 sublattices total
                        switch i
                            case 1
                                these_atoms = obj.l1_sublat1;
                                tagged_atoms = zeros(size(these_atoms,1),4);
                            case 2
                                these_atoms = obj.l1_sublat2;
                                tagged_atoms = zeros(size(these_atoms,1),4);
                            case 3
                                these_atoms = l2_sublat1_withwrap;
                                atom_storage = cell(5);
                            case 4
                                these_atoms = l2_sublat2_withwrap;
                                atom_storage = cell(5);
                        end
                        
                        for xx = 1:DIV_INCREMENT
                            for yy = 1:DIV_INCREMENT
                                thisxlims = [xcrit(xx),xcrit(xx+1)];
                                thisylims = [ycrit(yy),ycrit(yy+1)];
                                xtest = (these_atoms(:,1) <= thisxlims(1)) & (these_atoms(:,1) > thisxlims(2));
                                ytest = (these_atoms(:,2) <= thisylims(1)) & (these_atoms(:,2) > thisylims(2));
                                test = xtest & ytest;
                                if i < 3
                                    tagged_atoms(test,:) = [these_atoms(test,:),repmat([xx,yy],size(these_atoms(test,:),1),1)];
                                else
                                    atom_storage{xx,yy} = these_atoms(test,:);
                                end
                            end
                        end
                        switch i
                            case 1
                                l1_sublat1_tagged = tagged_atoms;
                            case 2
                                l1_sublat2_tagged = tagged_atoms;
                            case 3
                                l2_sublat1_storage = atom_storage;
                            case 4
                                l2_sublat2_storage = atom_storage;
                        end
                            
                    end
                    
                    % Now go to actually do the DSC computation.
                    for i = 1:2  % 2 sublattices
                        switch i
                            case 1
                                these_tagged = l1_sublat1_tagged;
                                this_compare = l2_sublat1_storage;
                                all_compare = l2_sublat1_withwrap;
                                this_DSC = zeros(size(these_tagged,1),2);
                            case 2
                                these_tagged = l1_sublat2_tagged;
                                this_compare = l2_sublat2_storage;
                                all_compare = l2_sublat2_withwrap;
                                this_DSC = zeros(size(these_tagged,1),2);
                                
                        end
                        for j = 1:size(these_tagged,1)
%                             j
                            tags = these_tagged(j,3:4);
                            to_compare = this_compare{tags(1),tags(2)};
                            d = sum((to_compare - these_tagged(j,1:2)).^2,2).^0.5;
                            [val,idx] = min(d);
                            if val > 1.2
                                to_compare = all_compare;
                                d = sum((to_compare - these_tagged(j,1:2)).^2,2).^0.5;
                                [~,idx] = min(d);
                            end
                            this_DSC(j,:) =  to_compare(idx,:) - these_tagged(j,1:2);
%                             num_to_compare = size(this_compare{tags(1),tags(2)},1);
%                             for k = 1:num_to_compare
%                                 this_compare{tags(1),tags(2)}(k,:);
                        end
                        switch i
                            case 1
                                DSC_sublat1 = this_DSC;
                            case 2
                                DSC_sublat2 = this_DSC;
                        end
                    end
                     
                case 'old'
                    %%% NPK 12/26/2019: changed this so that it is layer2 - layer1,
                    %%% which will make the DSC the vector pointing from layer1 to
                    %%% layer2 as intended. However, it still doesn't quite match
                    %%% the simulation fitting results and I'm not sure why.
                    for i = 1:numlat1
                        d = sum((l2_sublat1_withwrap - obj.l1_sublat1(i,:)).^2,2).^0.5;
                        [val,idx] = min(d);
                        DSC_sublat1(i,:) =  l2_sublat1_withwrap(idx,:) - obj.l1_sublat1(i,:);
%                         i
                    end
                    
                    numlat2 = size(obj.l1_sublat2,1);
                    DSC_sublat2 = zeros(size(obj.l1_sublat2));
                    for i = 1:numlat2
                        d = sum((l2_sublat2_withwrap - obj.l1_sublat2(i,:)).^2,2).^0.5;
                        [val,idx] = min(d);
                        DSC_sublat2(i,:) = l2_sublat2_withwrap(idx,:) - obj.l1_sublat2(i,:);
%                         i
                    end
            end
            
            obj.DSC_sublat1 = DSC_sublat1;
            obj.DSC_sublat2 = DSC_sublat2;
        end
        
        
        
        % Helper Function for assigning psuedostacking regions
        function filtered_lat = filter_DSC(obj,type,r)
            DSC_field = obj.getDSClat();
            lat = obj.getLayer1();
            t = 60/180*pi;
            rotmat = [cos(t) sin(t); -sin(t) cos(t)];
            v1 = [0;obj.a/sqrt(3)];
            v2 = rotmat*v1;
            if nargin < 3
                r = 0.3;
            end
            
            switch type
                case 'AA'
                    [tf] = isInCircle(DSC_field(:,1),DSC_field(:,2),0,0,r);
                    filtered_lat = lat(tf,:);
                case 'AB'
                    h(1,:) = v1;
                    h(2,:) = v2;
                    h(3,:) = v2-v1;
                    h(4,:) = -v1;
                    h(5,:) = -v2;
                    h(6,:) = v1-v2;
                    [tf1] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(1,1),h(1,2),r);
                    [tf2] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(2,1),h(2,2),r);
                    [tf3] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(3,1),h(3,2),r);
                    [tf4] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(4,1),h(4,2),r);
                    [tf5] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(5,1),h(5,2),r);
                    [tf6] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(6,1),h(6,2),r);
                    filtered_lat = lat(tf1|tf2|tf3|tf4|tf5|tf6,:);
                case 'SP'
                    h(1,:) = (v1+v2)/2;
                    h(2,:) = -(v1+v2)/2;
                    h(3,:) = (2*v1-v2)/2;
                    h(4,:) = -(2*v1-v2)/2;
                    h(5,:) = (2*v2-v1)/2;
                    h(6,:) = -(2*v2-v1)/2;
                    [tf1] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(1,1),h(1,2),r);
                    [tf2] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(2,1),h(2,2),r);
                    [tf3] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(3,1),h(3,2),r);
                    [tf4] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(4,1),h(4,2),r);
                    [tf5] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(5,1),h(5,2),r);
                    [tf6] = isInCircle(DSC_field(:,1),DSC_field(:,2),h(6,1),h(6,2),r);
                    filtered_lat = lat(tf1|tf2|tf3|tf4|tf5|tf6,:);
            end
            
            
        end
        
        
        
        % Use the above function to assign psuedostacking regions. Make
        % this a quantitative readout at some point.
        function [figh,perc_AA,perc_AB,perc_SP] = assignPsuedostacking(obj,figh)
            filtered_AAlat = obj.filter_DSC('AA');
            filtered_ABlat = obj.filter_DSC('AB');
            filtered_SPlat = obj.filter_DSC('SP');
            if nargin < 2
                figh = figure;
            else
                figure(figh);
            end
            hold on
            sh1 = scatter(filtered_AAlat(:,2),filtered_AAlat(:,1),12,'filled','r');
            sh2 = scatter(filtered_ABlat(:,2),filtered_ABlat(:,1),12,'filled','g');
            sh3 = scatter(filtered_SPlat(:,2),filtered_SPlat(:,1),12,'filled','b');
            legend([sh1,sh2,sh3],'AA','AB','SP');
            
            alldsc = obj.getDSClat();
            nDSC = size(alldsc,1);
            nAA = size(filtered_AAlat,1);
            nAB = size(filtered_ABlat,1);
            nSP = size(filtered_SPlat,1);
            perc_AA = nAA/nDSC*100;
            perc_AB = nAB/nDSC*100;
            perc_SP = nSP/nDSC*100;
        end
        
        % Visualize all of the DSC vector points that have been obtained
        function figh = plotDSCLattice(obj)
            DSC_lat = obj.getDSClat();
            figh = figure;
            scatter(DSC_lat(:,1),DSC_lat(:,2),10,'filled');
            axis equal
            hold on
            t = 60/180*pi;
            rotmat = [cos(t) sin(t); -sin(t) cos(t)];
            vec1 = [0;obj.a/sqrt(3)];
            vec2 = rotmat*vec1;
            mat = [vec1,vec2];
            quiver([0;0],[0;0],mat(1,:)',mat(2,:)',0,'r')
            xlabel('DSC vector x component');
            ylabel('DSC vector y component');
            title('DSC Reduced Lattice: Scatterplot of DSC Vectors');
        end
        
        
        function figh = plotDSCQuiver(obj,figh)
            if nargin == 1 || isempty(figh)
                figh = figure
                clf
                set(gcf,'color','w')
            else
                figure(figh);
                hold on
            end
            DSC_lat = obj.getDSClat();
            tblg_lat1 = obj.getLayer1();
            qh = quiver(tblg_lat1(:,2),tblg_lat1(:,1),DSC_lat(:,2),DSC_lat(:,1),0,'k');
            %             quiver(obj.l1_sublat2(:,1),obj.l1_sublat2(:,2),obj.DSC_sublat2(:,1),obj.DSC_sublat2(:,2));
            
            xlabel('y (Angstroms)');
            ylabel('x (Angstroms)');
            title(obj.name);
            line([0 0 obj.cellDimXY(2) obj.cellDimXY(2) 0],...
                [0 obj.cellDimXY(1) obj.cellDimXY(1) 0 0],...
                'linewidth',1,'color','k')
            axis equal
            set(gca,'ydir','reverse')
            box on
            xlim([0 obj.cellDimXY(2)] + [-2 2])
            ylim([0 obj.cellDimXY(1)] + [-2 2])
            set(gca,'position',[0.02 0.04 0.96 0.92])
            set(gcf, 'Position',  [0, 0, 500, 800])
        end
        
        % Visualize the DSC vector field. If not stored, then compute it first.
        function figh = plotDSCField(obj,figh)
            DSC_lat = obj.getDSClat();
            tblg_lat1 = obj.getLayer1();
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
            if nargin < 2
                figh = figure; 
            else
                figure(figh);
            end
            contourf(ybase,xbase,interpAmp',50,'LineStyle','None');
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
        
        
        function plotAveragedDSCLattice(obj,stepsize,averaging_radius)
            xsteps = 0:stepsize:obj.cellDimXY(1);
            ysteps = 0:stepsize:obj.cellDimXY(2);
            l1_coords = obj.getLayer1();
            DSCvectors = obj.getDSClat();
            averageDSC = zeros(numel(xsteps)*numel(ysteps),2);
            count = 1;
            for i = 1:numel(xsteps)
                for j = 1:numel(ysteps)
                    circle_center = [xsteps(i), ysteps(j)];
                    [tf] = isInCircle(l1_coords(:,1),l1_coords(:,2),circle_center(1),circle_center(2),averaging_radius);
                    averageDSC(count,:) = mean(DSCvectors(tf,:),1);
                    count = count + 1;
                end
            end
            
            figure
            scatter(averageDSC(:,1),averageDSC(:,2),10,'filled');
            axis equal
            hold on
            t = 60/180*pi;
            rotmat = [cos(t) sin(t); -sin(t) cos(t)];
            vec1 = [0;obj.a/sqrt(3)];
            vec2 = rotmat*vec1;
            mat = [vec1,vec2];
            quiver([0;0],[0;0],mat(1,:)',mat(2,:)',0,'r');
            axh = gca;
            plotFullDisplacementHexagons( axh );
            xlabel('DSC vector x component');
            ylabel('DSC vector y component');
            title('DSC Reduced Lattice: Scatterplot of DSC Vectors');
            
        end
        
        
        function blinking_predictions = predictBlinkingPatterns(obj,nan_handling_flag)
            if isempty(obj.DSC_sublat1)
                obj.computeDSCField();
            end
            
            NUM_DISKS = 12;
            DSC_lat = obj.getDSClat();
            scaling_constant = 1;
            blinking_predictions = zeros(size(DSC_lat,1),NUM_DISKS);
            for i = 1:size(DSC_lat,1)
                [ this_pred ] = trigFittingFunctions( DSC_lat(i,:), scaling_constant, nan_handling_flag );
                blinking_predictions(i,:) = this_pred;  % files all twelve into the columns of the matrix.
            end
        end
       
        
        % Actually make the blinking plots, greyscale.
        function makeBlinkingPlots(obj,subplot_flag,include_redundant_flag)
            if subplot_flag && include_redundant_flag
                warning('Cannot make subplots with all 12 redundant disks -- making separate figures instead.');
                subplot_flag = 0;
            end
            blinking_predictions = obj.predictBlinkingPatterns();
            l1coords = obj.getLayer1();  % because the DSC lattice directly corresponds to these
            % make the interpolation points
            npoints = 500;
            xbase = linspace(0,obj.cellDimXY(1),npoints);
            ybase = linspace(0,obj.cellDimXY(2),npoints+1);
            [xspace,yspace] = meshgrid(xbase,ybase);
            
            if include_redundant_flag
                disks_to_plot = 1:12;
            else
                disks_to_plot = [1:3,7:9];
            end
            
            if subplot_flag
                figure;
            end
            
            counter = 0;
            for i = disks_to_plot
                counter = counter + 1;
                if subplot_flag
                    subplot(2,3,counter);
                    set(gcf, 'Position',  [0, 0, 750, 600])
                else
                    figure;
                    set(gcf, 'Position',  [0, 0, 500, 800])
                end
                this_blink = blinking_predictions(:,i);
                myInterpBlink = scatteredInterpolant(l1coords(:,1),l1coords(:,2),this_blink,'linear','nearest');
                interpBlink = myInterpBlink(xspace,yspace);
                interpBlink(interpBlink > 1) = nan;
                contourf(ybase,xbase,interpBlink',50,'LineStyle','None');
                colormap(gray)
                axis equal
                shading interp
                colorbar
                title(sprintf(strcat(['Blinking prediction: disk ',num2str(i)])));
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
            end
        end
        
        
        
        % Plotting method for the tblg and offset blg classes specifically
        % figh is optional; if left blank, will make a new figure.
        % Otherwise, will do a hold on write onto the currently opened
        % figure.
        function figh = plot(obj,figh,new_convention)
            if nargin < 3
                new_convention = false;
            end
            if exist('figh','var')
                if ~isempty(figh)
                    figure(figh);
                    hold on;
                end
            else
                figh = figure;
                clf;
                hold on
                set(gcf,'color','w')
            end
            
            if new_convention
                scatter(obj.l1_sublat1(:,1),obj.l1_sublat1(:,2),'marker','.','sizedata',40,...
                    'markeredgecolor',[1 0 0],'markerfacecolor','none','linewidth',1);
                scatter(obj.l1_sublat2(:,1),obj.l1_sublat2(:,2),'marker','.','sizedata',40,...
                    'markeredgecolor',[0.7 0 0.3],'markerfacecolor','none','linewidth',1);
                
                scatter(obj.l2_sublat1(:,1),obj.l2_sublat1(:,2),'marker','o','sizedata',20,...
                    'markeredgecolor',[00 0.7 0.9],'markerfacecolor','none','linewidth',1);
                scatter(obj.l2_sublat2(:,1),obj.l2_sublat2(:,2),'marker','o','sizedata',20,...
                    'markeredgecolor',[00 0.4 1],'markerfacecolor','none','linewidth',1);
            else
                scatter(obj.l1_sublat1(:,2),obj.l1_sublat1(:,1),'marker','.','sizedata',40,...
                    'markeredgecolor',[1 0 0],'markerfacecolor','none','linewidth',1);
                scatter(obj.l1_sublat2(:,2),obj.l1_sublat2(:,1),'marker','.','sizedata',40,...
                    'markeredgecolor',[0.7 0 0.3],'markerfacecolor','none','linewidth',1);
                
                scatter(obj.l2_sublat1(:,2),obj.l2_sublat1(:,1),'marker','.','sizedata',30,...
                    'markeredgecolor',[00 0.7 0.9],'markerfacecolor','none','linewidth',1);
                scatter(obj.l2_sublat2(:,2),obj.l2_sublat2(:,1),'marker','.','sizedata',30,...
                    'markeredgecolor',[00 0.4 1],'markerfacecolor','none','linewidth',1);

%                 scatter(obj.l1_sublat1(:,2),obj.l1_sublat1(:,1),'marker','.','sizedata',40,...
%                     'markeredgecolor',[1 0 0],'markerfacecolor','none','linewidth',1);
%                 scatter(obj.l1_sublat2(:,2),obj.l1_sublat2(:,1),'marker','.','sizedata',40,...
%                     'markeredgecolor',[0.7 0 0.3],'markerfacecolor','none','linewidth',1);
                
%                 scatter(obj.l2_sublat1(:,2),obj.l2_sublat1(:,1),'marker','o','sizedata',20,...
%                     'markeredgecolor',[00 0.7 0.9],'markerfacecolor','none','linewidth',1);
%                 scatter(obj.l2_sublat2(:,2),obj.l2_sublat2(:,1),'marker','o','sizedata',20,...
%                     'markeredgecolor',[00 0.4 1],'markerfacecolor','none','linewidth',1);
            end
            
            legend('Layer 1, sublattice 1','Layer 1, sublattice 2',...
                'Layer 2, sublattice 1','Layer 2, sublattice 2');
            if new_convention
                xlabel('x (Angstroms)');
                ylabel('y (Angstroms)');
            else
                xlabel('y (Angstroms)');
                ylabel('x (Angstroms)');
            end
            title(obj.name);
            
            
            
            % line([0 uProj(2)],...
            %     [0 uProj(1)],...
            %     'linewidth',1,'color','g')
            % line([0 vProj(2)],...
            %     [0 vProj(1)],...
            %     'linewidth',1,'color','g')
            hold off
            axis equal
            if new_convention
                set(gca,'ydir','normal');
            else
                set(gca,'ydir','reverse')
            end
            
            box on
            
            if new_convention
                set(gca,'position',[0.02 0.04 0.96 0.92])
                set(gcf, 'Position',  [50, 50, 1000, 500])
                xlim([0 obj.cellDimXY(1)] + [-2 2])
                ylim([0 obj.cellDimXY(2)] + [-2 2])
                line([0 0 obj.cellDimXY(1) obj.cellDimXY(1) 0],...
                [0 obj.cellDimXY(2) obj.cellDimXY(2) 0 0],...
                'linewidth',1,'color','k')
            else
                set(gca,'position',[0.02 0.04 0.96 0.92])
                set(gcf, 'Position',  [0, 0, 500, 800])
                xlim([0 obj.cellDimXY(2)] + [-2 2])
                ylim([0 obj.cellDimXY(1)] + [-2 2])
                line([0 0 obj.cellDimXY(2) obj.cellDimXY(2) 0],...
                [0 obj.cellDimXY(1) obj.cellDimXY(1) 0 0],...
                'linewidth',1,'color','k')
            end
        end
        
        
        % Added by NPK on 10/30/2020
        % Will draw little circles instead of markers.
        function plotVolumetric(obj)
            l1 = obj.getLayer1();
            l2 = obj.getLayer2();
            l1 = horzcat(l1,zeros(size(l1,1),1));
            l2 = horzcat(l2,obj.c*ones(size(l2,1),1));
            
            % plot
            figh = figure;
            scatterSpherical( l1, 0.35, [0.4,0.4,0.4], figh );
            scatterSpherical( l2, 0.32, [0.722,0.451,0.20], figh );
%             camlight
%             lighting phong
            buffer = obj.cellDimXY(1)/5;
            pbaspect([1,1,1]);
            xlim([-buffer,obj.cellDimXY(1)+buffer]);
            ylim([-buffer,obj.cellDimXY(2)+buffer]);
            zlim([-buffer-obj.cellDimXY(1)/2,buffer+obj.cellDimXY(1)/2]);
        end
        
        
        
        % NPK 10/30/2020
        % Displaces all atoms by an interpolated displacement field. This
        % should be a phase-unwrapped, extended zone displacement field. 
        %
        % Goes equal and opposite on the two lattices. 3d array, first
        % component is the x displacement, second componet is the y
        % displacement.
        function splineDisplacement(obj,dfield,spacing)
             xdfield = dfield(:,:,1);
             ydfield = dfield(:,:,2);
             % convert from nm to Angstrom
             [xmesh,ymesh] = meshgrid(10*spacing*(spacing:size(xdfield,1)),10*spacing*(spacing:size(xdfield,2)));
             [newx_displacements_l1sublat1] = interp2(xmesh,ymesh,xdfield/2,obj.l1_sublat1(:,1),obj.l1_sublat1(:,2));
             [newx_displacements_l1sublat2] = interp2(xmesh,ymesh,xdfield/2,obj.l1_sublat2(:,1),obj.l1_sublat2(:,2));
             [newy_displacements_l1sublat1] = interp2(xmesh,ymesh,ydfield/2,obj.l1_sublat1(:,1),obj.l1_sublat1(:,2));
             [newy_displacements_l1sublat2] = interp2(xmesh,ymesh,ydfield/2,obj.l1_sublat2(:,1),obj.l1_sublat2(:,2));
             
             [newx_displacements_l2sublat1] = interp2(xmesh,ymesh,-xdfield/2,obj.l2_sublat1(:,1),obj.l2_sublat1(:,2));
             [newx_displacements_l2sublat2] = interp2(xmesh,ymesh,-xdfield/2,obj.l2_sublat2(:,1),obj.l2_sublat2(:,2));
             [newy_displacements_l2sublat1] = interp2(xmesh,ymesh,-ydfield/2,obj.l2_sublat1(:,1),obj.l2_sublat1(:,2));
             [newy_displacements_l2sublat2] = interp2(xmesh,ymesh,-ydfield/2,obj.l2_sublat2(:,1),obj.l2_sublat2(:,2));
             
             obj.l1_sublat1 = obj.l1_sublat1 + [newx_displacements_l1sublat1,newy_displacements_l1sublat1];
             obj.l1_sublat2 = obj.l1_sublat2 + [newx_displacements_l1sublat2,newy_displacements_l1sublat2];
             obj.l2_sublat1 = obj.l2_sublat1 + [newx_displacements_l2sublat1,newy_displacements_l2sublat1];
             obj.l2_sublat2 = obj.l2_sublat2 + [newx_displacements_l2sublat2,newy_displacements_l2sublat2];
        end
        
        
        
        function invertCoords(obj)
            obj.l1_sublat1 = -obj.l1_sublat1;
            obj.l1_sublat2 = -obj.l1_sublat2;
            obj.l2_sublat1 = -obj.l2_sublat1;
            obj.l2_sublat2 = -obj.l2_sublat2;
        end
    end
    
end

