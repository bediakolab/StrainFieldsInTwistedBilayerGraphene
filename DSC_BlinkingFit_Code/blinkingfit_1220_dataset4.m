% blinkingfit_1220_dataset4.m

filename = '4_50x50_ss=0p7nm_C2=40um_alpha=1mrad_spot7_100ms_CL=130 bin=4_scan90_60kV_NPK_Normalized.h5';
disk_filename = '3_50x50_ss=0p7nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPK_Normalized__integration_disks.mat';
md4 = FourDSTEM_Analysis_Engine(filename,disk_filename);
% [d,si,ei] = m4.partialLoad(10);
% size(d)
% si
% ei

md4.setBlinkingDisks();
md4.plotIntegrationDisks();
md4.integrateDisks();
md4.makeBlinkingPlots();
md4.fitBlinking();
