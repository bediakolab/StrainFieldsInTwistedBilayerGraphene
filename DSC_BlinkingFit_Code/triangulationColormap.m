function [ cmap ] = triangulationColormap()
% Colormap function for making the uniaxial heterostrain % faceted
% triangulation plots. 

len1 = 20;
len2 = 30;
len3 = 30;

% hval = 0.72;
hval = 0.8;
hval2 = 0.1;
% hval4 = 0.72;
hval4 = 0.72;
hval3 = 0.75;

valval = 0.4;
% hval = 0.05;
c4hsv = [hval4,0.2,1];
c3hsv = [hval3,0.7,1];
c2hsv = [hval,1,valval];
% c2hsvupper = [hval-1,1,valval];
c2hsvupper = [hval-1,1,valval];
c1hsv = [hval2,0,0];

% c1rgb = hsv2rgb(c1hsv);
% c2rgb = hsv2rgb(c2hsv);
% c3rgb = hsv2rgb(c3hsv);
c1rgb = c1hsv;
c2rgb = c2hsv;
c3rgb = c3hsv;
c4rgb = c4hsv;

mycmap = vertcat(c1rgb,c2hsvupper);
seg1 = [linspace(mycmap(1,1),mycmap(2,1),len1)',linspace(mycmap(1,2),mycmap(2,2),len1)',linspace(mycmap(1,3),mycmap(2,3),len1)'];
seg2 = [linspace(c2rgb(1),c3rgb(1),len2)',linspace(c2rgb(2),c3rgb(2),len2)',linspace(c2rgb(3),c3rgb(3),len2)'];
seg3 = [linspace(c3rgb(1),c4rgb(1),len3)',linspace(c3rgb(2),c4rgb(2),len3)',linspace(c3rgb(3),c4rgb(3),len3)'];


cmap = hsv2rgb(mod(vertcat(seg1,seg2,seg3),1.00001));

end

