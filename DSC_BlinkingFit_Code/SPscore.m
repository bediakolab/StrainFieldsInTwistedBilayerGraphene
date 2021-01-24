function [ scores ] = SPscore( reduced_displacements )
% Function that scores a displacement vector for degree of saddle point
% character; max is 1, min is zero. Only on the basis of the angle; no
% amplitude involved.
%
% Code adapted and simplified from getCustomDisplacementColor.m.

full_angles = atan2(reduced_displacements(:,2),reduced_displacements(:,1));
in1 = full_angles >= 0 & full_angles < pi/6;
in2 = full_angles >= pi/6 & full_angles < pi/3;
in3 = full_angles >= pi/3 & full_angles < pi/2;
in4 = full_angles >= pi/2 & full_angles < 2*pi/3;
in5 = full_angles >= 2*pi/3 & full_angles < 5*pi/6;
in6 = full_angles >= 5*pi/6 & full_angles <= pi;

AB_value = 0;
SP_value = 1;

H1 = full_angles/(pi/6)*AB_value + (1-full_angles/(pi/6))*SP_value;
H2 = (1-(full_angles-pi/6)/(pi/6))*AB_value + ((full_angles-pi/6)/(pi/6))*SP_value;
H3 = ((full_angles-pi/3)/(pi/6))*AB_value + (1-(full_angles-pi/3)/(pi/6))*SP_value;
H4 = (1-(full_angles-pi/2)/(pi/6))*AB_value + ((full_angles-pi/2)/(pi/6))*SP_value;
H5 = ((full_angles-2*pi/3)/(pi/6))*AB_value + (1-(full_angles-2*pi/3)/(pi/6))*SP_value;
H6 = (1-(full_angles-5*pi/6)/(pi/6))*AB_value + ((full_angles-5*pi/6)/(pi/6))*SP_value;

scores = H1.*in1 + H2.*in2 + H3.*in3 + H4.*in4 + H5.*in5 + H6.*in6;

end

