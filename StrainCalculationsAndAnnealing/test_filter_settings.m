% test_filter_settings.m
% 
% Script for trying to ascertain the extent of the impact of the filtering
% on the calculated strain parameters. Load the DS26 annealed, or whichever
% else.
%
% Nathanael Kazmierczak, 05/20/2020

moving_average_radii = 3:2:21;
n_filts = numel(moving_average_radii);

eyx_maxval_stor = zeros(n_filts,1);
eyx_minval_stor = zeros(n_filts,1);
eyx_avval_stor = zeros(n_filts,1);
eyx_percentiles_stor = zeros(101,n_filts);
gradientsum_percentiles_stor = zeros(101,n_filts);
eyx_values_stor = zeros(200^2,n_filts);
gradientsum_values_stor = zeros(200^2,n_filts);

for i = 1:n_filts
    i
    this_diameter = moving_average_radii(i);
    stringcell = {'median outlier','xdisp',[5,0.2];'median outlier','ydisp',[5,0.2];'moving average circle','xdisp',this_diameter;'moving average circle','ydisp',this_diameter};
    [ filterstruct ] = buildFilterStruct( stringcell );
    [eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(0,'flat',cmap,[50,0,0],[],[0,0],0,filterstruct);
    close all
    eyx_maxval = max(eyx_values(:));
    eyx_minval = min(eyx_values(:));
    eyx_avval = mean(eyx_values(:));
    eyx_maxval_stor(i) = eyx_maxval;
    eyx_minval_stor(i) = eyx_minval;
    eyx_avval_stor(i) = eyx_avval;
    eyx_percentiles_stor(:,i) = eyx_percentiles(:);
    eyx_values_stor(:,i) = eyx_values(:);
    gradientsum_values_stor(:,i) = gradientsum_values(:);
    gradientsum_percentiles_stor(:,i) = gradientsum_percentiles(:);
end

[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(0,'flat',cmap,[50,0,0],[],[0,0],0,filterstruct);

figure; histogram(eyx_values); title('eyx values for diameter = 21');
figure; histogram(gradientsum_values); title('gradientsum values for diameter = 21');
figure; plot(moving_average_radii,eyx_maxval_stor,'o-'); hold on;
plot(moving_average_radii,eyx_minval_stor,'o-');
plot(moving_average_radii,eyx_avval_stor,'o-');
legend('eyx maximum value','eyx minimum value','eyx average value');

figure
for i = 1:n_filts
    plot((0:100)',eyx_percentiles_stor(:,i),'-');
    hold on
end
xlabel('percentile');
ylabel('eyx strain value');
legend('3','5','7','9','11','13','15','17','19','21');

figure
for i = 1:n_filts
    plot((0:100)',gradientsum_percentiles_stor(:,i),'-');
    hold on
end
xlabel('percentile');
ylabel('gradientsum strain value');
legend('3','5','7','9','11','13','15','17','19','21');

f = figure;
dumpfig = figure;
for i = 1:n_filts
    figure(dumpfig);
    h = histogram(eyx_values_stor(:,i));
    centers = movmean(h.BinEdges,2);
    centers = centers(2:end);
    values = h.Values;
    figure(f)
    plot(centers(:),values(:),'-');
    hold on
end
ylabel('histogram bin counts');
xlabel('eyx strain value %');
legend('3','5','7','9','11','13','15','17','19','21');

f2 = figure;
for i = 1:n_filts
    figure(dumpfig);
    h = histogram(gradientsum_percentiles_stor(:,i));
    centers = movmean(h.BinEdges,2);
    centers = centers(2:end);
    values = h.Values;
    figure(f2);
    plot(centers(:),values(:),'-');
    hold on
end
ylabel('histogram bin counts');
xlabel('total gradient strain value %');
legend('3','5','7','9','11','13','15','17','19','21');

figure;
for i = 1:n_filts
    [f,xi] = ksdensity(eyx_values_stor(:,i));
    plot(xi(:),f(:),'-');
    hold on
end
ylabel('kernel density');
xlabel('eyx strain value %');
legend('3','5','7','9','11','13','15','17','19','21');

figure;
for i = 1:n_filts
    [f,xi] = ksdensity(gradientsum_percentiles_stor(:,i));
    plot(xi(:),f(:),'-');
    hold on
end
ylabel('kernel density');
xlabel('total gradient strain value %');
legend('3','5','7','9','11','13','15','17','19','21');
