% get_cropped_DS1S1_newps.m
%
% Nathanael Kazmierczak, 01/23/2021

fp = '/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed/06022020_chiralfixed_DS1_postannealing.mat';
load(fp);

trim_value = 0;
stringcell = {'median outlier','both displacements',[3,1]};
pre_filterstruct = buildFilterStruct(stringcell);
TGVid = 'DS1';
calibration_method = 'SP';
divide_by_two_for_intralayer = 'false';
sign_convention = [1,1];
trim_tear_value = [];
SPdirection = 1;
m4.computeStrainMaps3(trim_value,pre_filterstruct,TGVid,calibration_method,divide_by_two_for_intralayer,sign_convention,trim_tear_value,SPdirection);

% Crop the displacements 
m4.dfield_filtered_for_strain = m4.dfield_filtered_for_strain(10:100,100:190,:);
use_annealed = 2; % New setting, makes other things unused because references operations by computeStrainMaps3.
croprange = [];  % unused
filterstruct = [];
%     r = 0.3;
%     r = 'modified wigner-seitz';
r = 'intermediate';
make_plots = true;
[stacking_struct,~,~,pseudostacking_mat] = m4.getPseudostackingAssignments(make_plots,r,filterstruct,use_annealed,croprange);



ps_storage = [stacking_struct.AAinnerperc,stacking_struct.AAouterperc,...
    stacking_struct.SPperc,stacking_struct.ABperc,stacking_struct.ABSPtransitionperc];
