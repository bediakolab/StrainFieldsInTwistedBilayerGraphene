% interactive_blinking_script.m
%
% Script for showing how the predicted blinking pattern varies as a
% function of the displacement vector.
%
% Nathanael Kazmierczak, 02/29/2020

[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
% t = 60/180*pi;
% rotmat = [cos(t) sin(t); -sin(t) cos(t)];
% v1 = [0,hexagon_lattice_constant/sqrt(3)];
% v2 = v1*rotmat';
% vl1 = vectorLine(v2-v1,v1);
% vl2 = vectorLine(v1,v2);
% vl3 = vectorLine(v2,-v1);
% vl4 = vectorLine(v2-v1,-v1);
% vl5 = vectorLine(v1,-v2);
% vl6 = vectorLine(-v2,v1);


figh = figure;
s1 = subplot(2,2,1);
hold on
axis equal
plotLine( figh, v1, v2 );
plotLine( figh, v2, v2-v1 );
plotLine( figh, -v1, v2-v1 );
plotLine( figh, -v1, -v2 );
plotLine( figh, v1-v2, -v2 );
plotLine( figh, v1-v2, v1 );

% [ pred_vals ] = trigFittingFunctions( DSC_guess, scaling_constant, bound_handling_flag );

plotLine(figh, 2*v2, 3*v2 - v1);
plotLine(figh, 3*v2 - 2*v1, 3*v2 - v1);
plotLine(figh, 3*v2 - 2*v1, 2*v2 - 2*v1);
plotLine(figh, 2*v2 - 2*v1, 2*v2 - 2*v1);
plotLine(figh, 2*v2 - 2*v1, v2-v1);

plotLine(figh, v1, 2*v1);
plotLine(figh, 2*v1 + v2, 2*v1);
plotLine(figh, 2*v1 + v2, 2*v2 + v1);
plotLine(figh, 2*v2, 2*v2 + v1);
plotLine(figh, 2*v2, v2);

plotLine(figh, v1 - 3*v1, 2*v1 - 3*v1);
plotLine(figh, 2*v1 + v2 - 3*v1, 2*v1 - 3*v1);
plotLine(figh, 2*v1 + v2 - 3*v1, 2*v2 + v1 - 3*v1);
plotLine(figh, 2*v2 - 3*v1, 2*v2 + v1 - 3*v1);
plotLine(figh, 2*v2 - 3*v1, v2 - 3*v1);
plotLine(figh, -2*v1, -3*v1 + v2);

plotLine(figh, v1 - 3*v2, 2*v1 - 3*v2);
plotLine(figh, 2*v1 + v2 - 3*v2, 2*v1 - 3*v2);
plotLine(figh, 2*v1 + v2 - 3*v2, 2*v2 + v1 - 3*v2);
plotLine(figh, 2*v2 - 3*v2, 2*v2 + v1 - 3*v2);
plotLine(figh, 2*v2 - 3*v2, v2 - 3*v2);
plotLine(figh, -2*v2, -3*v2 + v1);

plotLine(figh, v1 - 3*v2 + v1 + v2, 2*v1 - 3*v2 + v1 + v2);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2, 2*v1 - 3*v2 + v1 + v2);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2, 2*v2 + v1 - 3*v2 + v1 + v2);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2, 2*v2 + v1 - 3*v2 + v1 + v2);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2, v2 - 3*v2 + v1 + v2);
plotLine(figh, -2*v2 + v1 + v2, -3*v2 + v1 + v1 + v2);

plotLine(figh, v1 - 3*v2 + v1 + v2 - 3*v1, 2*v1 - 3*v2 + v1 + v2 - 3*v1);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2 - 3*v1, 2*v1 - 3*v2 + v1 + v2 - 3*v1);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2 - 3*v1, 2*v2 + v1 - 3*v2 + v1 + v2 - 3*v1);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2 - 3*v1, 2*v2 + v1 - 3*v2 + v1 + v2 - 3*v1);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2 - 3*v1, v2 - 3*v2 + v1 + v2 - 3*v1);
plotLine(figh, -2*v2 + v1 + v2 - 3*v1, -3*v2 + v1 + v1 + v2 - 3*v1);

plotLine( figh, v1, v2, 'r' );
plotLine( figh, v2, v2-0.5*v1, 'r' );
plotLine( figh, v2-0.5*v1, -v2+0.5*v1, 'r' );
plotLine( figh, -v2+0.5*v1, -v2+v1, 'r' );
plotLine( figh, v1-v2, v1, 'r' );
hold on;

clear('impointCallback');
color1 = [1,0,0];
color2 = [0,0,0];
h = impoint();
h2 = impoint();
h.setColor(color1);
h2.setColor(color2);
s2 = subplot(2,2,2);
s3 = subplot(2,2,3);
s4 = subplot(2,2,4);
axis equal
fcn = @(pos) impointCallback(pos,h2,s2,s1,s3,s4,color1,1);
% fcn(h.getPosition());
fcn2 = @(pos) impointCallback(pos,h,s2,s1,s3,s4,color2,2);
fcn(h.getPosition());
fcn2(h2.getPosition());
id = addNewPositionCallback(h,fcn);
id2 = addNewPositionCallback(h2,fcn2);





