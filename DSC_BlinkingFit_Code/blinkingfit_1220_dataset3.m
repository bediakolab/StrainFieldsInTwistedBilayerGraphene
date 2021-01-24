% blinkingfit_1220_dataset3.m
%
% First script with the new analysis class.

filename = '3_50x50_ss=0p7nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPK_Normalized.h5';
addpath(genpath(pwd));

m4 = FourDSTEM_Analysis_Engine(filename);
% [d,si,ei] = m4.partialLoad(10);
% size(d)
% si
% ei

m4.setBlinkingDisks();
m4.integrateDisks();
m4.makeBlinkingPlots();
m4.fitBlinking();
