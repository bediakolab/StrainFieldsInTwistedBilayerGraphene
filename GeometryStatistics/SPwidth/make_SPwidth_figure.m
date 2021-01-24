% make_SPwidth_figure.m
%
% Script to make an error bar plot. This script uses the 'average of the
% ranges' calculation technique, which was the second one developed.
% 
% A version with both standard error and 95% confidence intervals is
% computed (since the current 95% confidence intervals are astronomical).

twist_angle = [1.37
1.23
1.23
1.03
0.63
0.26
1.17
0.75
0.633
0.6632
0.16];

twist_std = [0.0315
0.0284
0.024
0.0333
0.0189
0.0248
0.0444
0.0299
0.0239
0.0229
0.03];

SPrange_average = [1.967
2.2542
2.2384
2.4027
3.4522
3.2795
2.2601
3.0965
3.2416
3.1827
3.6179];

stderrs = [0.0411
0.0405
0.0456
0.0552
0.165
0.3491
0.0565
0.1064
0.3505
0.2059
0.468];

conf95 = [0.0822
0.081
0.0912
0.1104
0.33
0.6982
0.113
0.2128
0.701
0.4118
0.936];

figure;
errorbar(twist_angle,SPrange_average,stderrs,stderrs,twist_std,twist_std,'o');
xlabel('Twist angle (degrees)');
ylabel('SP width (nm)');
title('Using standard error for uncertainty quantification.');

figure;
errorbar(twist_angle,SPrange_average,conf95,conf95,twist_std,twist_std,'o');
xlabel('Twist angle (degrees)');
ylabel('SP width (nm)');
title('Using 95% confidence interval for uncertainty quantification.');


