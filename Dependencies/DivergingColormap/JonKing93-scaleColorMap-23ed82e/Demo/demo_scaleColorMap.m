% This is a demo of scaleColorMap

%% Figure 1: Scale temperature anomalies (a diverging dataset) to a
% blue-white-red colormap

% Load a blue-white-red colormap. This scales from blue (for the smallest
% values), to white (at the center of the colormap), to red (for maximum
% values).
m = load('demo_colormap.mat');
cmap = m.cmap;

% Load some demo temperature anomaly data. (Specifically, temperature
% anomaly data over the Americas for December, July, and October relative
% to the annual mean temperature.)
m = load('demo_T_anomaly_data.mat');
Tanom = m.Tanom;

% Plot the data
figure;
h = pcolor(Tanom(:,:,1));
set(h,'linestyle','none');
colorbar
set(gca,'xtick',[],'ytick',[]);
title('October Temperature Anomaly (C)');

% Using scaleColorMap adjusts the colormap and color limits so that the
% white values at the center of the colormap are matched to the anomaly
% divergence point of zero.
%
% Furthermore, this data is from boreal winter. So, the
% negative anomalies are much more extreme than the positive anomalies.
% Accordingly, scaleColorMap adjusts the colorbar to use more intense blue
% values, while less intense red values.
scaleColorMap( cmap, 0 );



%% Figure 2: Scale and use the same colormap for multiple plots.

% It's often useful to compare several plots while using the same colorbar.
% Let's do that here...

% Plot the data
figure

ax1 = subplot(1,3,1);
h = pcolor(Tanom(:,:,1));
set(h,'linestyle','none');
colorbar
set(gca,'xtick',[],'ytick',[]);
title('December Anomaly (C)');

ax2 = subplot(1,3,2);
h = pcolor(Tanom(:,:,2));
set(h,'linestyle','none');
colorbar
set(gca,'xtick',[],'ytick',[]);
title('July Anomaly (C)');

ax3 = subplot(1,3,3);
h = pcolor(Tanom(:,:,3));
set(h,'linestyle','none');
colorbar
set(gca,'xtick',[],'ytick',[]);
title('October Anomaly (C)');

% Get the set of axes handles
allAxes = [ax1, ax2, ax3];

% Scale the data on all the plots to the same colorbar
scaleColorMap( cmap, 0, allAxes );
