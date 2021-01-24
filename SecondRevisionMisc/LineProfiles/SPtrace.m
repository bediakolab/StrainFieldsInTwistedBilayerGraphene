% SPtrace.m
%
% 01/19/2021 Nathanael Kazmierczak
% Use twist angle 0.26 degrees: DS26S1
%
% It doesn't appear that I saved any interlayer rotational plots in the
% main folder. 

addpath('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
addpath('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2');
load('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
dataset_cell = '05312020_chiralfixed_noprefactor_DS26_postannealing.mat';

if false
load('01192021_SPtracedata.mat');
figure; 
xdata = xdata-11.75;
plot(xdata,ydata);
ylabel('Local reconstruction rotation (degrees)');
xlabel('Distance from SP center (nm)');
savedata = [xdata',ydata'];
xlswrite('01192021_SProtationprofile_data.xlsx',savedata);
end


if false
load(dataset_cell);
% very light outlier filtering to start off
stringcell = {'median outlier','both displacements',[3,0.5]};
pre_filterstruct = buildFilterStruct(stringcell);
TGVID = 'DS26';  % versus Colin's original parameter settings.
rot_calibration_id = 'SP';
tensor_rotate_SPid = 1;
colormaptype = {cmap,fire};  % option to have one colormap for principal strain.
divide_by_two_for_intralayer = false;
sign_convention = [+1,-1];
% sign_convention = [0,0];  % This would be the correct comparison to Muller: divide by 2 but do not remove Moire
trim_tear_value = 0;
m4.computeStrainMaps3(30,pre_filterstruct,TGVID,rot_calibration_id,divide_by_two_for_intralayer,sign_convention,trim_tear_value,tensor_rotate_SPid);

saveplotflag = false;
overlay_registration = false;
m4.plotStrainMaps3(saveplotflag,overlay_registration,colormaptype);



m4.dfield_filtered_for_strain;
figure;
imagesc(m4.dfield_filtered_for_strain(:,:,1)); set(gca,'yDir','normal'); axis equal; title('xdisp');
figure
imagesc(m4.dfield_filtered_for_strain(:,:,2)); set(gca,'yDir','normal'); axis equal; title('ydisp');
end

if false

load('ImprofileSPdisplacementData.mat');
xdisp_offset = 15;
ydisp_offset = 15;
xdisp_xdata = xdisp_xdata/2 - xdisp_offset;  % The original scale wasn't set to nm for the x and y filtered displacement graphs, but rather pixels. 
ydisp_xdata = ydisp_xdata/2 - ydisp_offset;
xdisp_ydata = xdisp_ydata*10;
ydisp_ydata = ydisp_ydata*10;


figure;
plot(xdisp_xdata,xdisp_ydata);
hold on
plot(ydisp_xdata,ydisp_ydata);
xlabel('Distance from SP region (nm)');
ylabel('Displacement (Angstrom)');
end
% load(dataset_cell{1});
% m4.make

%%%Actually it turns out this is really annoying because would have to get
%%%the two data series to line up. Let's just pick a vertical saddle point
%%%instead. 


load('SecondImprofileSPDataBetter.mat')
tensiledisp_offset = 14;
sheardisp_offset = 14;
tensiledisp_xdata = tensiledisp_x/2 - tensiledisp_offset;  % The original scale wasn't set to nm for the x and y filtered displacement graphs, but rather pixels. 
sheardisp_xdata = sheardisp_x/2 - sheardisp_offset;
tensiledisp_ydata = tensiledisp_y*10;
sheardisp_ydata = sheardisp_y*10;

verticaloffset = sheardisp_ydata(1);
sheardisp_ydata = -(sheardisp_ydata - verticaloffset);


figure;
plot(tensiledisp_xdata,tensiledisp_ydata);
hold on
plot(sheardisp_xdata,sheardisp_ydata);
xlabel('Distance from SP region (nm)');
ylabel('Displacement (Angstrom)');
