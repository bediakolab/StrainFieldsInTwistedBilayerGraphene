% hotspot_probability_calculation.m


addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset4MaddieFit');
load('objectdataDS4toaddSP.mat');
m4.filename = 'TEAM1Day1Dataset4DiffractionSI.h5';
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset4MaddieFit';
m4.saved_plots_foldername = 'HotspotProbability';

%% Visualizations, histograms
SEC1 = false;
if SEC1
    hs_indices = [60,102];
    nhs_indices = [60,103];
    % m4.plotSingleDP(hs_indices(1),hs_indices(2),make_histograms_flag,remove_hBN,remove_Graphene,remove_Beamstop,correct_DP);
    m4.plotSingleDP(hs_indices(1),hs_indices(2),1,0,0,0,0);
    m4.plotSingleDP(hs_indices(1),hs_indices(2),0,1,0,1,0);
    [data] = m4.singleLoad(hs_indices(1),hs_indices(2));
    figure
    edges = -0.25:0.5:10.25;
    histogram(data(:),edges);
end
% all that poisson MaxLik fit to ROI stuff

%% The proper simulation
% (1) Assume all pixels are blinked fully on in the AA region. Calculate
% Poisson probability of getting lower.

m4.getNumBlinkingDiskPixels();
niter = 1024;  % power of 2

DSC_xstorage = zeros(niter,1);
DSC_ystorage = zeros(niter,1);
for i = 1:niter
    i
    integration_vals = zeros(1,12);
    true_displacement = [0,0];
    expected_trig_blinking = trigFittingFunctions( true_displacement, m4.trig_prefactors, 0 );
    for q = 1:12
%         trig_prefactor = m4.trig_prefactors(q);
%         numpixels = m4.disk_pixel_nums(q);
        numpixels = 4;
        r = poissrnd(expected_trig_blinking(q),1,numpixels);
        integration_vals(q) = mean(r);
    end
    fitfun = @(DSC_guess) integration_vals - trigFittingFunctions( DSC_guess, m4.trig_prefactors, 0 );
    ub = [1.24,1.43];
    lb = [-1.24,0];
    trigfun_flag = 1;
    options = optimoptions('lsqnonlin');
    options.Display = 'off';
    [ DSC_results ] = multistartDiskOptimization( fitfun, lb, ub, options, trigfun_flag );
    DSC_xstorage(i) = DSC_results(1);
    DSC_ystorage(i) = DSC_results(2);
end
displacement_field = cat(3,reshape(DSC_xstorage,[sqrt(niter),sqrt(niter)]),reshape(DSC_ystorage,[sqrt(niter),sqrt(niter)]));
[r,c,h] = size(displacement_field);

figure;
x = displacement_field(:,:,1);
y = displacement_field(:,:,2);
scatter(x(:),y(:),'filled');
axh = gca;
plotFullDisplacementHexagons(axh);



% Copied out of the custom color plot generation
[ RGB_color_stack ] = getCustomDisplacementColor( displacement_field );
f = figure;
set(f,'Position',[200,200,800,700]);
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);
axes(ax1);
imagesc(1:r,1:c,RGB_color_stack);
set(ax1,'ydir','normal');
axis square
xlabel('Real space x (nm)','FontSize',12);
ylabel('Real space y (nm)','FontSize',12);
title('Colorized AB, SP, and AA regions','FontSize',14);
hold on
set(ax1,'Position',[0.05 0.1 0.7 0.75]);
set(ax1,'FontSize',12);

grid_density = 0.01;
[ RGB_color_stack_legend ] = getTriangleColorLegend(grid_density);
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


