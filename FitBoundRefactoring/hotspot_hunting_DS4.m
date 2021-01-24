% hotspot_hunting_DS4.m

if ~exist('m4','var')
    load('objectdataDS4toaddSP.mat');
end
% m4.makeCustomDisplacementColorPlot
indices = [59,101] + 1;  % because the plots always start from zero nm displacement
% note that the pixels both above and below are also hotspots.
hotspot_displacement = permute(m4.DSC_fit_storage(indices(1),indices(2),:),[1,3,2])
hotspot_integration_values = permute(m4.disk_averages(indices(1),indices(2),:),[3,1,2])
%  It looks like signal averaging in displacement space will make this go
%  away. Bin 2 is most crude, but could also weight by 5-10% from the
%  nearest datasets and that might be enough.
% hotspot_integration_values = permute(m4.disk_averages(indices(1):indices(1)+1,indices(2):indices(2)+1,:),[3,1,2])
% hotspot_integration_values = mean(mean(hotspot_integration_values,3),2)

nonhs_indices = [60,103];
nonhs_displacement = permute(m4.DSC_fit_storage(nonhs_indices(1),nonhs_indices(2),:),[1,3,2])
nonhs_integration_values = permute(m4.disk_averages(nonhs_indices(1),nonhs_indices(2),:),[3,1,2])

weight_vector_local = ones(1,12);
objfcn_hs = @(DSCvector) (hotspot_integration_values' - trigFittingFunctions(DSCvector,m4.trig_prefactors,0)).*weight_vector_local;
objfcn_nhs = @(DSCvector) (nonhs_integration_values' - trigFittingFunctions(DSCvector,m4.trig_prefactors,0)).*weight_vector_local;


xbase = -1.5:0.01:1.5;
ybase = 0:0.01:1.5;
[xspace,yspace] = meshgrid(xbase,ybase);
RMSR_hs = zeros(size(xspace));
RMSR_nhs = zeros(size(xspace));
residuals_hs = zeros([size(xspace),12]);
residuals_nhs = zeros([size(xspace),12]);
for i = 1:numel(ybase)
    i
    for j = 1:numel(xbase)
        this_DSC = [xspace(i,j),yspace(i,j)];
        result_hs = objfcn_hs(this_DSC);
        result_nhs = objfcn_nhs(this_DSC);
        RMSR_hs(i,j) = rms(result_hs);
        RMSR_nhs(i,j) = rms(result_nhs);
        residuals_hs(i,j,:) = permute(result_hs,[1,3,2]);
        residuals_nhs(i,j,:) = permute(result_nhs,[1,3,2]);
    end
end

% RMSR_hs(isnan(RMSR_hs)) = 0;
% RMSR_nhs(isnan(RMSR_nhs)) = 0;

figure
subplot(2,1,1);
contourf(xbase,ybase,RMSR_hs,30);
title('Hotspot RMSR plot (Dataset 4, indices [60,102])','FontSize',12)
axis equal
colormap jet
colorbar
caxis([0.5,3]);
subplot(2,1,2);
contourf(xbase,ybase,RMSR_nhs,30);
title('Non-hotspot RMSR plot (Dataset 4, indices [60,103])','FontSize',12)
colormap jet
axis equal 
colorbar
caxis([0.5,3]);

for i = 1:12
    figure
    subplot(2,1,1);
    contourf(xbase,ybase,residuals_hs(:,:,i),10);
    title(sprintf('Hotspot disk %d residuals plot (Dataset 4, indices [60,102])',i),'FontSize',12)
    axis equal
    colormap jet
    colorbar
    caxis([0,5]);
    
    subplot(2,1,2);
    contourf(xbase,ybase,residuals_nhs(:,:,i),10);
    title(sprintf('Non-hotspot disk %d residuals plot (Dataset 4, indices [60,103])',i),'FontSize',12)
    colormap jet
    axis equal
    colorbar
    caxis([0,5]);
end

rms_inner_disk_residuals_hs = rms(residuals_hs(:,:,1:6),3);
rms_outer_disk_residuals_hs = rms(residuals_hs(:,:,7:12),3);
rms_inner_disk_residuals_nhs = rms(residuals_nhs(:,:,1:6),3);
rms_outer_disk_residuals_nhs = rms(residuals_nhs(:,:,7:12),3);

figure
subplot(2,1,1);
contourf(xbase,ybase,rms_inner_disk_residuals_hs,20);
title(sprintf('Hotspot inner ring RMSR plot (Dataset 4, indices [60,102])'),'FontSize',12)
axis equal
colormap jet
colorbar
caxis([0,4]);
subplot(2,1,2);
contourf(xbase,ybase,rms_inner_disk_residuals_nhs,20);
title(sprintf('Non-hotspot inner ring RMSR plot (Dataset 4, indices [60,103])'),'FontSize',12)
colormap jet
axis equal
colorbar
caxis([0,4]);

figure
subplot(2,1,1);
contourf(xbase,ybase,rms_outer_disk_residuals_hs,20);
title(sprintf('Hotspot outer ring RMSR plot (Dataset 4, indices [60,102])'),'FontSize',12)
axis equal
colormap jet
colorbar
caxis([0,4]);
subplot(2,1,2);
contourf(xbase,ybase,rms_outer_disk_residuals_nhs,20);
title(sprintf('Non-hotspot outer ring RMSR plot (Dataset 4, indices [60,103])'),'FontSize',12)
colormap jet
axis equal
colorbar
caxis([0,4]);
