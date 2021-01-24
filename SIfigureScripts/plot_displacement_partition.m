% plot_displacement_partition.m
%
% Nathanael Kazmierczak, 07/22/2020

SPcolor =  [0.9373, 0.3039, 0];
AAcolor = [0,0,0];
ABcolor = [1,1,1];

[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
v3 = v1-v2;
% patch vertices are from the bottom left going CCW
mp1 = (v1+v3)/2;
mp2 = (v1+v2)/2;
pv = [v3(1),0;
      v3(1),v3(2)/2;
      v3;
      (mp1+v3)/2;
      (mp1+v1)/2;
      v1;
      (mp2+v1)/2;
      (mp2+v2)/2;
      v2;
      v2(1),v2(2)/2;
      v2(1),0];

figure; 
c1 = [0,0;
      pv(1:2,:)];
patch(c1(:,1),c1(:,2),SPcolor);
hold on
c2 = [0,0;
     pv(2:4,:)];
patch(c2(:,1),c2(:,2),ABcolor);
c3 = [0,0;
     pv(4:5,:)];
patch(c3(:,1),c3(:,2),SPcolor);
c4 = [0,0;
     pv(5:7,:)];
patch(c4(:,1),c4(:,2),ABcolor);
c5 = [0,0;
     pv(7:8,:)];
patch(c5(:,1),c5(:,2),SPcolor);
c6 = [0,0;
     pv(8:10,:)];
patch(c6(:,1),c6(:,2),ABcolor);
c7 = [0,0;
     pv(10:11,:)];
patch(c7(:,1),c7(:,2),SPcolor);
axis equal

r = 0.71;
N = 1000;
t = linspace(0,pi,N);
x = r*cos(t);
y = r*sin(t);
patch(x,y,AAcolor);
xlim([-1.4,1.4]);
xlabel('x displacement (Angstrom)');
ylabel('y displacement (Angstrom)');
title('Pseudostacking displacement partition');
% optional addition "v2"
if true
    plotFullDisplacementHexagons(gca);
    xlim([-1.4,1.4]);
    ylim([-0.3,1.7]);
end
