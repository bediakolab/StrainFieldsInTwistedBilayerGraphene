function plotStrain(strain_struct,titlestring,whichlayer,threshold_filter)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here
dim = size(strain_struct.exx_map);
xaxis = 1:dim(2);
yaxis = 1:dim(1);

UPPER_THRESHOLD = 0.005;
LOWER_THRESHOLD = -0.005;

if threshold_filter
    strain_struct.exx_map(strain_struct.exx_map > UPPER_THRESHOLD | strain_struct.exx_map < LOWER_THRESHOLD) = nan;
    strain_struct.eyy_map(strain_struct.eyy_map > UPPER_THRESHOLD | strain_struct.eyy_map < LOWER_THRESHOLD) = nan;
    strain_struct.exy_map(strain_struct.exy_map > UPPER_THRESHOLD | strain_struct.exy_map < LOWER_THRESHOLD) = nan;
    strain_struct.theta_map(strain_struct.theta_map > UPPER_THRESHOLD | strain_struct.theta_map < LOWER_THRESHOLD) = nan;
end

subplot(2,4,1 + (whichlayer-1)*4);
pcolor(xaxis,yaxis,strain_struct.exx_map)
colormap(jet);
colorbar;
xlabel('Dimension 4 of datacube');
ylabel('Dimension 3 of datacube');
legend('exx');
title(titlestring);
shading flat
hold on;

subplot(2,4,2 + (whichlayer-1)*4);
pcolor(xaxis,yaxis,strain_struct.eyy_map)
colormap(jet);
colorbar;
xlabel('Dimension 4 of datacube');
ylabel('Dimension 3 of datacube');
legend('eyy')
shading flat
title(titlestring);

subplot(2,4,3 + (whichlayer-1)*4);
pcolor(xaxis,yaxis,strain_struct.exy_map)
colormap(jet);
colorbar;
xlabel('Dimension 4 of datacube');
ylabel('Dimension 3 of datacube');
legend('exy')
shading flat
title(titlestring);

subplot(2,4,4 + (whichlayer-1)*4);
pcolor(xaxis,yaxis,strain_struct.theta_map)
colormap(jet);
colorbar;
xlabel('Dimension 4 of datacube');
ylabel('Dimension 3 of datacube');
legend('theta')
shading flat
title(titlestring);

end

