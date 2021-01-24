% trblg_linesearch_driver2.m
%
% Script for performing detailed multistart simulations to ascertain the
% functional form of the blinking.
%
% This script probes multiple lines in x, y, and diagonal directions.
%
% Nathanael Kazmierczak, Dec 2019


input('HAAAAALLLLLTTTTT! Make sure you are not about to overwrite data.');
% savename = '12272019_trblg_linesearches2_data2.mat';
savename = 'test.mat';
load('12202019multislice_disk_centers.mat');
moire_cell_to_match = 1.2;
v2_x = 1.2305;
hexagon_lattice_constant = 2.461;
v1_y = hexagon_lattice_constant/sqrt(3);
v2_y = 0.7104;
trblg_raster = cell(100,3);

% According to my calculation, these are the boundaries of the DSC hexagon
% in the horizontal direction.
linenum = 3;
linepoints = 100;
lines = cell(1,linenum);
lines{1} = horzcat(zeros(linepoints,1),linspace(-v1_y,v1_y,linepoints)');

x0 = 0;
y0 = v1_y;
m = (v2_y - y0)/(v2_x);
lineform_2_y = @(x) m.*(x-x0) + y0;
xvals = linspace(0,v2_x,linepoints)';
lines{2} = horzcat(xvals,lineform_2_y(xvals));


xvals3 = linspace(-v2_x/2,v2_x/2,linepoints)';
x03 = 0;
y03 = 0;
m03 = lineform_2_y(v2_x/2) / (v2_x/2);
lineform_3_y = @(x) m03.*(x-x03) + y03;
lines{3} = horzcat(xvals3,lineform_3_y(xvals3));

% 12 when using '12202019multislice_disk_centers.mat' file
disks_used = size(disk_centers,1);
averages = zeros(linepoints,disks_used,linenum);
mytimer = tic;
for k = 1:linenum
    thisline = lines{k};
    for i = 1:linepoints
        fprintf('Beginning multislice iteration %d out of %d for line %d of %d\n',i,linepoints,k,linenum);
        fprintf('\tElapsed time is %f seconds\n',toc(mytimer));
        thisDSCvector = thisline(i,:);
        trblg_raster{i,k} = TranslatedBilayerGraphene(thisDSCvector,moire_cell_to_match);
        trblg_raster{i,k}.simulate();
        trblg_raster{i,k}.setAveragingCircles(disk_centers,radius);
        thisaverages = trblg_raster{i,k}.getAverageDiskIntensities();
        averages(i,:,k) = thisaverages';
        save(savename);
    end
end

averages = averages - min(averages);
averages = averages./repmat(max(averages),size(averages,1),1);

trblg_raster{1,1}.plotAveragingCircles()

% figure;
% hold on
% scatter(linesearch_yvalue_points,averages(:,1));
% xlabel('DSC y component');
% ylabel('Blinking intensity, disk 1');
% myfun = @(y) 0.85*abs(cos(y*(pi/2)/v2_x));
% funvals = myfun(linesearch_yvalue_points);
% plot(linesearch_yvalue_points,funvals,'-r');
% myfun2 = @(y) cos(y*(pi/2)/v2_x).^2;
% funvals2 = myfun2(linesearch_yvalue_points);
% plot(linesearch_yvalue_points,funvals2,'-k');
% legend('Multislice Points','0.85*abs(cos)','cos^2');

figure
subplot(2,3,1);
hold on;
scatter(lines{1}(:,2),averages(:,1,1));
title('Disk 1');
xlabel('Line 1 y component');
subplot(2,3,2);
scatter(lines{1}(:,2),averages(:,2,1));
title('Disk 2');
xlabel('Line 1 y component');
subplot(2,3,3);
scatter(lines{1}(:,2),averages(:,3,1));
title('Disk 3');
xlabel('Line 1 y component');
subplot(2,3,4);
scatter(lines{1}(:,2),averages(:,7,1));
title('Disk 7');
xlabel('Line 1 y component');
subplot(2,3,5);
scatter(lines{1}(:,2),averages(:,8,1));
title('Disk 8');
xlabel('Line 1 y component');
subplot(2,3,6);
scatter(lines{1}(:,2),averages(:,9,1));
title('Disk 9');
xlabel('Line 1 y component');
% mtit('Line 1: DSC Lattice Vector 1');
% legend('Disk 1','Disk 3','Disk 7','Disk 8');


figure
hold on
subplot(2,3,1);
scatter(lines{2}(:,1),averages(:,1,2));
title('Disk 1')
xlabel('Line 2 x component');
subplot(2,3,2);
scatter(lines{2}(:,1),averages(:,2,2));
title('Disk 2')
xlabel('Line 2 x component');
subplot(2,3,3);
scatter(lines{2}(:,1),averages(:,3,2));
title('Disk 3')
xlabel('Line 2 x component');
subplot(2,3,4);
scatter(lines{2}(:,1),averages(:,7,2));
title('Disk 7')
xlabel('Line 2 x component');
subplot(2,3,5);
scatter(lines{2}(:,1),averages(:,8,2));
title('Disk 8')
xlabel('Line 2 x component');
subplot(2,3,6);
scatter(lines{2}(:,1),averages(:,9,2));
title('Disk 9')
xlabel('Line 2 x component');
% mtit('Line 2: Connect lattice vectors across saddle point');
% legend('Disk 1','Disk 2','Disk 3','Disk 7','Disk 8','Disk 9');

figure
hold on
subplot(2,3,1);
scatter(lines{3}(:,1),averages(:,1,3));
title('Disk 1')
xlabel('Line 3 x component');
subplot(2,3,2);
scatter(lines{3}(:,1),averages(:,2,3));
title('Disk 2')
xlabel('Line 3 x component');
subplot(2,3,3);
scatter(lines{3}(:,1),averages(:,3,3));
title('Disk 3')
xlabel('Line 3 x component');
subplot(2,3,4);
scatter(lines{3}(:,1),averages(:,7,3));
title('Disk 7')
xlabel('Line 3 x component');
subplot(2,3,5);
scatter(lines{3}(:,1),averages(:,8,3));
title('Disk 8')
xlabel('Line 3 x component');
subplot(2,3,6);
scatter(lines{3}(:,1),averages(:,9,3));
title('Disk 9')
xlabel('Line 3 x component');
% mtit('Line 3: DSC Lattice Vector 2');

movienum = 100;
DP1 = trblg_raster{1,1}.DP;
stack3D = zeros([size(DP1),movienum]);
for i = 1:movienum
    stack3D(:,:,i) = trblg_raster{i,3}.DP;
end
ftitle = '12272019line3movie_SP_betweenDSC1andDSC2_toSP';
writeMovie01NPK(stack3D,ftitle);

