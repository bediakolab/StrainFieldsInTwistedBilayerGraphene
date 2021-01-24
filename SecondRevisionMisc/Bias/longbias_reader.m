% longbias_reader.m
%
% For figuring out the effect of beam averaging.
% Distinguish from the beamwidth_testing_script2.m because here we are
% doing indicatrix circles in displacement space.
%
% Nathanael Kazmierczak, 03/23/2020

%% Set points for the simulation
% [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
% % grid_density = [51,51];
% grid_density = [10,10];
% [ raster_points ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant);
% raster_points(raster_points(:,2) < 0,:) = [];


addpath('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/SensitivityBiasData');
load('long_bias_simulation.mat');
f = figure;

%% Obtain the averaged blinking patterns.
radius = 0.2;
hold on
% viscircles(raster_points,repmat(radius,size(raster_points,1),1),'Color',[0.9,0.9,1]);
circlecolor = [0.8,0,0];

finalcolor = [0.7,0.8,0.6]; 
scattercolor = [0.9,0.6,0.9];
scatter(raster_points(:,1),raster_points(:,2),'filled','MarkerFaceColor',scattercolor,'MarkerEdgeColor',0.8*scattercolor);
color = 'k';
line_thickness = 2;
plotHalfDisplacementHexagon( gca, line_thickness, color );
% N = 1000;
% N = 1000;
% xbase = linspace(-2,2,N);
% ybase = linspace(-2,2,N);
% [xspace,yspace] = meshgrid(xbase,ybase);
% scaling_constants = ones(1,12);
% bound_handling_flag = 1; % not eliminating alleged nans
% fitted_displacement_field = zeros(size(raster_points,1),2);
% options.Display = 'off';
% lb = [-1.24,0];
% trigfun_flag = 1;
% ub = [1.24,1.43];
% for i = 1:size(raster_points,1)
%     fprintf('Analyzing raster point %d of %d.\n',i,size(raster_points,1));
%     this_center = raster_points(i,:);
%     [tf] = isInCircle(xspace,yspace,this_center(1),this_center(2),radius);
%     xvals = xspace(tf);
%     yvals = yspace(tf);
%     pred_val_stor = zeros(size(xvals,1),12);
%     for q = 1:size(xvals,1)
%         this_val = [xvals(q), yvals(q)];
%         pred_val_stor(q,:) = trigFittingFunctions( this_val, scaling_constants, 1 );
%     end
%     mean_pred_vals = mean(pred_val_stor,1);
%     
%     % Now fit this mean blinking pattern
%     objfcn = @(DSC_guess) mean_pred_vals - trigFittingFunctions( DSC_guess, scaling_constants, 0 );
%     [fitted_displacement_field(i,:),convergence_locations] = extendedMultistartDiskOptimization( objfcn, lb, ub, options, trigfun_flag, false );
% end

%% Plot biasing results
figure(f)
hold on
% finalcolor = [0.9,0.4,0];

scatter(fitted_displacement_field(:,1),fitted_displacement_field(:,2),'filled','MarkerFaceColor',finalcolor,'MarkerEdgeColor',finalcolor*0.8);
differences = fitted_displacement_field-raster_points;
magnitudes = (differences(:,1).^2 + differences(:,2).^2).^0.5;
quiver(raster_points(:,1),raster_points(:,2),differences(:,1),differences(:,2),0,'k')%,'Linewidth',1);
xlim([-1.3,1.3]);
ylim([-0.05,1.5]);
set(gca,'YDir','normal');
xlabel('x displacement (Angstrom)');
ylabel('y displacement (Angstrom)');

viscircles([0,0.71],radius,'Color',circlecolor);
scatter(0,0.71,'filled','MarkerFaceColor',circlecolor,'MarkerEdgeColor',0.8*circlecolor);

% 
% figure
% si = scatteredInterpolant(raster_points(:,1),raster_points(:,2),magnitudes,'nearest','none');
% xb = linspace(-v2(1)+0.01,v2(1)-0.01,grid_density(1));
% yb = linspace(0+0.01,v1(2)-0.01,grid_density(2));
% [xq,yq] = meshgrid(xb,yb);
% vq = si(xq,yq);
% imagesc(xb,yb,vq);
% colormap jet
% colorbar
% title(sprintf('Biasing magnitude for displacement width %.2f',radius));
% hold on
% plotFullDisplacementHexagons(gca);
% xlim([-1.3,1.3]);
% ylim([-0.05,1.5]);
% set(gca,'YDir','normal');

