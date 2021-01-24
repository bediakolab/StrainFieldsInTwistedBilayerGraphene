% makeAAradiusplot.m
%
% This uses the Gaussian circle data.
%
% Values are copied from the workbook
% 04252020GeometryStatisticsTabulationNPK_AACirclefit.xlsx
% and 05232020_TwistAngleAnalysisSummary.xlsx.
%
% Nathanael Kazmierczak, 07/05/2020

twist_angle_data = [1.37
1.23
1.23
1.19
1.03
0.84
0.75
0.63
0.29
0.32
0.16];

twist_angle_std = [0.0315
0.0284
0.024
0.0383
0.0333
0.0547
0.0299
0.0189
0.0066
0.0121
0.03];

AAradius_mean = [2.2677
2.3681
2.4003
2.3668
2.5655
2.6403
2.6812
2.6326
2.7475
2.8417
2.897];

AAradius_95conf = [0.0231
0.0309
0.0232
0.0374
0.0351
0.0816
0.0891
0.1062
0.1099
0.3041
0.176];

figure
errorbar(twist_angle_data,AAradius_mean,AAradius_95conf,AAradius_95conf,twist_angle_std,twist_angle_std,'o');
xlabel('Moire twist angle (degrees)');
ylabel('AA radius (nm)');
title('Gaussian Circle AA estimation');

