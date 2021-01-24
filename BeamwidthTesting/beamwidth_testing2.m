% beamwidth_testing2.m
%
% Script for ascertaining the effects of blinking pattern averaging.
% Nathanael Kazmierczak

moire_angle_deg = 1.4;
recon_angle = 0.5;
recon_distance = 50;
use_old_reconfun = 0;
tblg = TwistedBilayerGraphene(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun);
new_alg_flag = 1;
tblg.computeDSCField(new_alg_flag);

%% Compute all blinking patterns
nan_handling_flag = 1;  % will compute extended zone
blinking_predictions = tblg.predictBlinkingPatterns(nan_handling_flag);

%% Get averages
stepsize = 10;
averaging_radius = 100;
xsteps = 0:stepsize:tblg.cellDimXY(1);
ysteps = 0:stepsize:tblg.cellDimXY(2);
l1_coords = tblg.getLayer1();
averaged_blinking_patterns = zeros(numel(xsteps)*numel(ysteps),12);
averaged_displacement_fit = zeros(numel(xsteps)*numel(ysteps),2);
count = 1;
for i = 1:numel(xsteps)
    i
    for j = 1:numel(ysteps)
        circle_center = [xsteps(i), ysteps(j)];
        [tf] = isInCircle(l1_coords(:,1),l1_coords(:,2),circle_center(1),circle_center(2),averaging_radius);
        averaged_blinking_patterns(count,:) = mean(blinking_predictions(tf,:),1);
        count = count + 1;
    end
end

%% Perform preprocessing
% mins = min(min(averaged_blinking_patterns));
% maxes = repmat(max(max(averaged_blinking_patterns)),[numel(xsteps)*numel(ysteps),12]);
% disk_average_storage = averaged_blinking_patterns./maxes;
disk_average_storage = averaged_blinking_patterns;

%% Fit
trig_prefactor = 0.95;
nan_handle_flag = 0;
count = 1;
weight_vector = ones(1,12);
fitted_displacement_field = zeros(numel(xsteps)*numel(ysteps),2);
options.Display = 'off';
lb = [-1.24,0];
trigfun_flag = 1;
ub = [1.24,1.43];  % really a/sqrt(3) and a/2, I believe.
for i = 1:numel(xsteps)
    i
    for j = 1:numel(ysteps)
        this_disk_average_storage = disk_average_storage(count,:);
        objfcn = @(DSCvector) (this_disk_average_storage - trigFittingFunctions(DSCvector,trig_prefactor,nan_handle_flag)).*weight_vector;
        fitted_displacement_field(count,:) = multistartDiskOptimization( objfcn, lb, ub, options, trigfun_flag );
        count = count + 1;
    end
end

%% Make plots
figure
scatter(fitted_displacement_field(:,1),fitted_displacement_field(:,2),10,'filled');
axis equal
hold on
[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
mat = [v1',v2'];
quiver([0;0],[0;0],mat(1,:)',mat(2,:)',0,'r');
axh = gca;
plotFullDisplacementHexagons( axh );
xlabel('DSC vector x component');
ylabel('DSC vector y component');
title('DSC Reduced Lattice: Scatterplot of DSC Vectors');


