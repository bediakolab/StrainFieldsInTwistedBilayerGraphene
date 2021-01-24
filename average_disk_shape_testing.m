% average_disk_shape_testing.m
%
% Driver script for verifying the algorithm in getAveragedDiskShape.m

addpath('C:\Users\Bediako Lab\Desktop\4DSTEM data');
addpath('C:\Users\Bediako Lab\Desktop\4DSTEM data\Maddiedata07292019');
addpath('C:\Users\Bediako Lab\Desktop\PythonWorkspace\py4DSTEM\py4DSTEM\process\utils');

addpath('/Volumes/SANDISK/BediakoLabLexar');
addpath('/Volumes/SANDISK/BediakoLabLexar/4DSTEMdata_updating');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM');

format long;

info = h5info('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5');
info.Groups.Groups(1).Groups(1).Groups.Datasets(1).Dataspace;  % looks like it is 70 x 8

start = [1,1,1,1];
count = [512,512,70,8];  % chunks are in the fourth dimension
% filename = '3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um.h5';
filename = '3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5';
this_data = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
firstDP = this_data(:,:,1,1);

% make into 3D array
DPs = merge4thdim( this_data );

%% Registration and indexing (standard, from before)
[ graphene1_lattice_storage, graphene2_lattice_storage, hk_graphene1_storage, hk_graphene2_storage, disk_all_storage ] = registerDatacube( DPs );

%% Extract the disks to be averaged.
% via external knowledge
niter = size(disk_all_storage,1);
disk_positions_cell = cell(niter,1);
for i = 1:niter
    graphene_intensities = 100*ones(size(disk_all_storage{i},1),1);
    origin = mean(disk_all_storage{i}(:,2:3),1,'omitnan');
    [first_graphene_disks,second_graphene_disks] = sortGraphenePeaks(disk_all_storage{i},graphene_intensities,origin);
    disk_positions_cell{i} = first_graphene_disks;
end

%% test the desired function




index_choice = [1,0];  % arbitrarily
window_width = 20;
% Get averaged probe
subpixel = 'False';
average_probe = getAveragedDiskShape(DPs,disk_positions_cell,hk_graphene1_storage,index_choice,window_width,subpixel);


