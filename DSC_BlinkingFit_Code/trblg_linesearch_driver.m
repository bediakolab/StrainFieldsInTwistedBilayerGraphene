% trblg_linesearch_driver.m
%
% Script for performing detailed multistart simulations to ascertain the
% functional form of the blinking.
% Nathanael Kazmierczak, Dec 2019


input('HAAAAALLLLLTTTTT! Make sure you are not about to overwrite data.');
savename = '12262019_trblg_linesearch.mat';
load('12202019multislice_disk_centers.mat');
moire_cell_to_match = 1.2;
hexagon_horizontal_boundary = 1.2305;
clear trblg_raster;

% According to my calculation, these are the boundaries of the DSC hexagon
% in the horizontal direction.
linesearch_yvalue_points = linspace(-1.2305,1.2305,100)';

% 12 when using '12202019multislice_disk_centers.mat' file
disks_used = size(disk_centers,1);
averages = zeros(size(linesearch_yvalue_points,1),disks_used);
mytimer = tic;
for i = 1:size(linesearch_yvalue_points,1)
    fprintf('Beginning multislice iteration %d out of %d\n',i,size(linesearch_yvalue_points,1));
    fprintf('\tElapsed time is %f seconds\n',toc(mytimer));
    thisDSCvector = [0,linesearch_yvalue_points(i)];
    trblg_raster(i) = TranslatedBilayerGraphene(thisDSCvector,moire_cell_to_match);
    trblg_raster(i).simulate();
    trblg_raster(i).setAveragingCircles(disk_centers,radius);
    thisaverages = trblg_raster(i).getAverageDiskIntensities();
    averages(i,:) = thisaverages';
    save(savename);
end

averages = averages - min(averages);
averages = averages./repmat(max(averages),size(averages,1),1);

trblg_raster(1).plotAveragingCircles()

figure;
hold on
scatter(linesearch_yvalue_points,averages(:,1));
xlabel('DSC y component');
ylabel('Blinking intensity, disk 1');
myfun = @(y) 0.85*abs(cos(y*(pi/2)/hexagon_horizontal_boundary));
funvals = myfun(linesearch_yvalue_points);
plot(linesearch_yvalue_points,funvals,'-r');
myfun2 = @(y) cos(y*(pi/2)/hexagon_horizontal_boundary).^2;
funvals2 = myfun2(linesearch_yvalue_points);
plot(linesearch_yvalue_points,funvals2,'-k');
legend('Multislice Points','0.85*abs(cos)','cos^2');

