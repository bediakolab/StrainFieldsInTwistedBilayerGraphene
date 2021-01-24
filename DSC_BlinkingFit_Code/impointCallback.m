function impointCallback(position,otherh,plothandle,plothandle2,ph3,ph4,color,myid)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

persistent ph ch
if isempty(ph)
    ph = zeros(1,2);
end
if isempty(ch)
    ch = [];
end
RADIUS = 0.1;

[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
axes(plothandle);
bound_handling_flag = 1;
scaling_constant = 1;
[ pred_vals ] = trigFittingFunctions( position, scaling_constant, bound_handling_flag );

positions = [v1;v2;v2-v1;-v1;-v2;v1-v2;...
             2*v1-v2;v1+v2;2*v2-v1;v2-2*v1;-v1-v2;v1-2*v2];
n = 500;
xbase = linspace(-3,3,n);
ybase = linspace(-3,3,n);
[xspace,yspace] = meshgrid(xbase,ybase);
r = 0.2;
plotmat_randimine = zeros(size(xspace));

for i = 1:12
    [tf] = isInCircle(xspace,yspace,positions(i,1),positions(i,2),r); 
    plotmat_randimine(tf) = pred_vals(i);
end

colormap(fire);
imagesc(plotmat_randimine);
shading flat
colorbar
caxis([0,1]);
% xlim([0,500]);
% ylim([0,500]);
axis square

basis = permn(-3:3,2);
basis(basis(:,1) == 0 & basis(:,2) == 0,:) = [];
% basis = [2 2
%          2 1 
%          2 0 
%          2 -1
%          2 -2
%          1 2
%          1 1
%          1 0 
%          1 -1
%          1 -2
%          0 2
%          0 1
%          
%          0 -1
%          0 -2
%          -1 2
%          -1 1
%          -1 0
%          -1 -1
%          -1 -2
%          -2 2
%          -2 1
%          -2 0 
%          -2 -1
%          -2 -2];

w1 = v1 + v2;
w2 = 2*v2 - v1;
W = [w1',w2'];
genpoints = W*basis';
genpoints_plus = genpoints + repmat(position',1,size(genpoints,2));
genpoints_minus = genpoints - repmat(position',1,size(genpoints,2));
genpoints_all = [genpoints_plus,genpoints_minus,-position'];

axes(plothandle2);
hold on
if ~(ph(myid) == 0)
    delete(ph(myid));
end
delete(ch);
ph(myid) = plot(genpoints_all(1,:),genpoints_all(2,:),'o','MarkerEdgeColor',color);
hold on
% ch = viscircles(position,RADIUS);
xlim([-4,4]);
ylim([-4,4]);

% Plot a random one of the red points on ph3.
axes(ph4);
randimine = randi(size(genpoints_all,2));
randpos = genpoints_all(:,randimine)';
plotmat_randimine = zeros(size(xspace));
[ pred_vals_randimine ] = trigFittingFunctions( randpos, scaling_constant, bound_handling_flag );
for i = 1:12
    [tf] = isInCircle(xspace,yspace,positions(i,1),positions(i,2),r); 
    plotmat_randimine(tf) = pred_vals_randimine(i);
end
imagesc(plotmat_randimine);
shading flat
colorbar
caxis([0,1]);
% xlim([0,500]);
% ylim([0,500]);
axis square
title(sprintf('DP2 index: %d',randimine));


axes(ph3);

otherposition = otherh.getPosition();
% d = DisplacementEquivalenceClass(position);
% do = DisplacementEquivalenceClass(otherposition);
% displacement = d - do;
% axes(plothandle2);
% hold on
% quiver(position(1),position(2),displacement(1),displacement(2),0);

end

