function [ dict_function, v ] = getLatticePlaneNormalVectors(symflag)
% This function outputs a function handle that is to be used as a
% dictionary. Input the disk number you want (according to the standard
% ordering), and it pops out the vector.
%
% Nathanael Kazmierczak, Dec 2019

hexagon_lattice_constant = 2.461;

% t = 60/180*pi;

if symflag
    syms a t
    DSCv1 = [0; a/sqrt(3)];
    assume(a>0);
    DSCv1 = simplify(DSCv1/norm(DSCv1));
    rotmat = @(in) subs([cos(t) sin(t); -sin(t) cos(t)],t,in);  % rotates clockwise.
else
    DSCv1 = [0,hexagon_lattice_constant/sqrt(3)]';
    % Normalize everything from the start here.
    DSCv1 = DSCv1/norm(DSCv1);
    rotmat = @(t) [cos(t) sin(t); -sin(t) cos(t)];  % rotates clockwise.
end


% set 12 lattice planes
v = zeros(2,12);
v(:,1) = rotmat(pi/6)*DSCv1;
v(:,2) = rotmat(3*pi/6)*DSCv1;
v(:,3) = rotmat(5*pi/6)*DSCv1;
v(:,4) = -rotmat(7*pi/6)*DSCv1;
v(:,5) = -rotmat(9*pi/6)*DSCv1;
v(:,6) = -rotmat(11*pi/6)*DSCv1;
v(:,7) = DSCv1;
v(:,8) = rotmat(pi/3)*DSCv1;
v(:,9) = rotmat(2*pi/3)*DSCv1;
v(:,10) = -rotmat(3*pi/3)*DSCv1;
v(:,11) = -rotmat(4*pi/3)*DSCv1;
v(:,12) = -rotmat(5*pi/3)*DSCv1;

v = rotmat(pi/2)*v;   % I screwed it up the first time.

v = v';  % Actually the way I did it above was the wrong sign for the way 
% that I am thinking about the diffraction spots. So flip them here.
% v(:,2) = -v(:,2);  % No, that's not right

dict_function = @(disk_number) v(disk_number,:);




end

