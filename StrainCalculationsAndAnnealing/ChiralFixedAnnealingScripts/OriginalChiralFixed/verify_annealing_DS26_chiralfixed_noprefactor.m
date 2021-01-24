% verify_annealing_DS26_2.m
%
% DS26 has already been annealed, but the strain values were not saved
% along with it. Here we will treat it with the same filters as all of the
% other 05/23/2020 annealings, as well as give it the manual triangulation.
%
% Nathanael Kazmierczak, 05/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/ChiralFixed');
load('05312020_chiralfixed_noprefactor_DS26_forannealing.mat');  
load('demo_colormap.mat');
m4.saved_plots_folderpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/ChiralFixed';
m4.saved_plots_foldername = 'DS26_noprefactor_newsettings_06032020';

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
sign_convention = [+1,-1];
trim_edge_pixel_num = [40,0,0];
trim_tear_pixel_num = [];
overlay_registration = false;
divide_by_two_for_intralayer = true;
SP_num_for_xaxis = 1;
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention,trim_tear_pixel_num,...
                overlay_registration,divide_by_two_for_intralayer,SP_num_for_xaxis);
            
overlay_registration = true;          
[eyx_percentiles,gradientsum_percentiles,eyx_values,gradientsum_values] = ...
                m4.makeStrainMapsFromBlinking2(saveplotflag,cmap,trim_edge_pixel_num,filterstruct,sign_convention,trim_tear_pixel_num,...
                overlay_registration,divide_by_two_for_intralayer,SP_num_for_xaxis);
            
        
        