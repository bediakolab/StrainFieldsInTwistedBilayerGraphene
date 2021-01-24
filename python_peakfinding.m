% python_peakfinding.m
%
% You must start Matlab out of the Anaconda prompt for this to work.
% Simply type "matlab" in the Anaconda prompt, and it would appear that
% Matlab will be able to start to access all of the libraries that Python
% can access.
%
% Function that shows how we can 

py.importlib.import_module('utils');

% myans_nosubpixel =  py.utils.get_maxima_2D(peaks,pyargs('subpixel','False'))
myans_subpixel = py.utils.get_maxima_2D(peaks,pyargs('subpixel',logical(1)))
myans_nosubpixel = py.utils.get_maxima_2D(peaks,pyargs('subpixel',logical(0)))
% Now these come in as numpy nD arrays
peak1 = myans_subpixel{1}
peak2 = myans_subpixel{2}
peak3 = myans_subpixel{3}
% So convert to lists
peak1_list = peak1.tolist();
peak2_list = peak2.tolist();
peak3_list = peak3.tolist();

figure;
contourf(peaks,50);
hold on;
% so the first peak that comes out is first array dimension, thus 
scatter(peak2_list{1}+1,peak1_list{1}+1,'r','filled');
scatter(peak2_list{2}+1,peak1_list{2}+1,'r','filled');
scatter(peak2_list{3}+1,peak1_list{3}+1,'r','filled');
