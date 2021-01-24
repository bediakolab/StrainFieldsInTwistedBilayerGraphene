% verify_annealing_DS26_2.m
%
% DS26 has already been annealed, but the strain values were not saved
% along with it. Here we will treat it with the same filters as all of the
% other 05/23/2020 annealings, as well as give it the manual triangulation.
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets');
load('dataset26_ver2_deterministic_annealed_04272020.mat');  
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults';
m4.saved_plots_foldername = 'DS26_05232020';

% Filteramp and filterrange should be populated by the loaded .mat file.
filteramp = 0.3;
filterrange = 9;
m4.setEmitterOrientations();
m4.annealDisplacementField2(filteramp,filterrange,1);

% Special to DS 26, need to perform a manual triangulation to get the
% correct adjacency.
TRIANGULATE = false;  % Don't want to have to go through the manual process every single time.
if TRIANGULATE
    m4.triangulateAA(1,true);
    m4.setFacetTwistAngles();
    [average_angle,std_angle,max_angle,min_angle,n_triangles_used] = m4.plotFacetedTwistAngles([0,inf])
end

method = 'outside input';
% AA_type_string = 'AA Gaussian ellipse';
% window_size = 80;
% window_increment = 10;
AA_type_string = [];
window_size = [];
window_increment = [];
crop_val_pixels = [];
angle_input = 0.2609;
m4.computeMoireAngleEstimates(method,AA_type_string,window_size,window_increment,crop_val_pixels,angle_input);

%% This works decently well we know

% As a good first approximation
saveplotflag = 1;
m4.makeCustomDisplacementColorPlot([],[],filteramp,filterrange,1,0,0,saveplotflag);
filterstruct = buildFilterStructPreset('DS26 annealed displacement #2');
sign_convention = [-1,-1];
trim_edge_pixel_num = [40,0,0];
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention);
