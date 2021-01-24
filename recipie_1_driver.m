% recipie_1_driver.m


%% Get data

addpath /Volumes/Lexar
info = h5info('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5');
h5disp('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data');

% Attempt to read in exactly three diffraction patterns
start = [1,1,1,1];
count = [512,512,20,1];
data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
firstDP = data(:,:,1);

%% Define the window masks through ginput
hexagonal_window_masks = getHexagonalWindow(firstDP);

%% Invoke analysis on first DP for testing.
kernel_radius = 5;
sigma = 0;
min_disk_distance = 8;
num_to_find_per_window = 2;
make_figures = 1;
lower_threshold = 0.5;
upper_threshold = 2;
[graphene_locations, graphene_intensities] = braggDiskRegistrationRecipie_1(firstDP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold);

%% Loop over other DPs to see if the window and xcorr2 algorithm holds up.
make_figures = 1;
for i = 2:size(data,3)
    i
    thisDP = data(:,:,i);
    [graphene_locations, graphene_intensities] = braggDiskRegistrationRecipie_1(thisDP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold);
end

