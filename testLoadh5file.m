% testLoadh5file.m
%
% Tests the loading of the .h5 files that can be saved out of py4DSTEM

% datasave is 08052019test_diffraction_patterns_with_beamstop.mat

addpath /Volumes/Lexar
info = h5info('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5');
h5disp('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data');

% Attempt to read in exactly three diffraction patterns
start = [1,1,1,1];
count = [512,512,3,1];
data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
% imshow(uint8(data(:,:,1)));
% figure
% contourf(data(:,:,1),200);

figure;
pcolor(data(:,:,1)); shading flat;

% in order to do imshow, the data must be stepped down to uint8 even though
% it has full uint16 precision.

% Data saved in 

% Note that info.Groups.Groups(1) is data, while (2) is log and (3) is
% metadata.

% info.Groups.Groups(1).Groups(1) is /4DSTEM_experiment/data/datacubes.
% info.Groups.Groups(1).Groups(2) is /4DSTEM_experiment/data/diffractionslices
% info.Groups.Groups(1).Groups(3) is /4DSTEM_experiment/data/pointlistarrays
% info.Groups.Groups(1).Groups(4) is /4DSTEM_experiment/data/pointlists
% info.Groups.Groups(1).Groups(5) is /4DSTEM_experiment/data/realslices

% When displaying
% info.Groups.Groups(1).Groups(1).Groups.Datasets(1).Dataspace,
% the dimensions are given as Size: [512 512 70 8]. This means that the
% first two dimensions are dimensions in reciprocal space (perhaps y and
% x?), the third is y (because that's the dimension that I left uncropped)
% and the fourth is x (because that's the dimension that I cropped.
