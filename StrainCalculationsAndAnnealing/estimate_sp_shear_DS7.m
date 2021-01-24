% estimate_sp_shear_DS7.m
%
% First attempt at getting the shear results through a method comparable to
% that of Muller's 2013 PNAS paper. 
%
% Nathanael Kazmierczak

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset7RadialResidualsInterp');
load('03232020MagneticAnnealingExtendedZoneReconstruction.mat');
clear('figh1');
close all;


% annealing results.
fa = figure
quiver(rgrid(:),cgrid(:),dfieldx(:),dfieldy(:),0);
hold on
quiver(fixed_emitters_indices(:,2),fixed_emitters_indices(:,1),fixed_emitters_directions(:,1),fixed_emitters_directions(:,2),0,'r');
set(gca,'YDir','normal');

% custom plot
m4.makeDisplacementMapsFromBlinking(0,[],0,'flat')

[fc,axc] = m4.makeCustomDisplacementColorPlot();

%% set boundaries and directions.
SET = false;
figure(4)
if SET

disp('Please select the left and right endpoints of the SP region.');
coords = ginput(2);
r = input('Please enter the distance in each direction the fitting region may go: ');
end

r = 13;

Linep = coords(2,:) - coords(1,:);
parallel_angle = atan(Linep(2)/Linep(1))
perp_angle = pi/2 + parallel_angle
hold on
inc = [r*cos(perp_angle),r*sin(perp_angle)];
bcoords1 = coords(1,:) + inc;
bcoords2 = coords(1,:) - inc;
bcoords3 = coords(2,:) + inc;
bcoords4 = coords(2,:) - inc;

line([bcoords1(1),bcoords2(1)],[bcoords1(2),bcoords2(2)],'Color','r');
line([bcoords2(1),bcoords4(1)],[bcoords2(2),bcoords4(2)],'Color','r');
line([bcoords4(1),bcoords3(1)],[bcoords4(2),bcoords3(2)],'Color','r');
line([bcoords3(1),bcoords1(1)],[bcoords3(2),bcoords1(2)],'Color','r');

figure(fc)
axes(axc)
line([bcoords1(1),bcoords2(1)],[bcoords1(2),bcoords2(2)],'Color','r');
line([bcoords2(1),bcoords4(1)],[bcoords2(2),bcoords4(2)],'Color','r');
line([bcoords4(1),bcoords3(1)],[bcoords4(2),bcoords3(2)],'Color','r');
line([bcoords3(1),bcoords1(1)],[bcoords3(2),bcoords1(2)],'Color','r');

% Linecut function to get x coords. base_t is a linear interpolation
% between bounding box edges 2 and 4.
linecut = @(base_t,N) [linspace(0,r*2,N)'.*cos(perp_angle),linspace(0,r*2,N)'.*sin(perp_angle)] + repmat((1-base_t)*bcoords2 + base_t*bcoords4,N,1);

%% Get a linecut average over the soliton
dfieldx_relative = dfieldx;  % because we are referencing this against the v1
dfieldy_relative = dfieldy; % change this momentarily
Ncuts = 50;
NinCut = 100;
linecutdisps = zeros(NinCut,2,Ncuts);
basecoords = linspace(0,1,Ncuts);
% xbase = 0.5*(1:size(dfieldx,1));
% ybase = 0.5*(1:size(dfieldy,1));
[xspace,yspace] = meshgrid(xbase,ybase);
for i = 1:Ncuts
    this_linecut = linecut(basecoords(i),NinCut);
%     dispxq = interp2(dfieldx,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
%     dispyq = interp2(dfieldy,xspace,yspace,this_linecut(:,1),this_linecut(:,2));
    dispxq = interp2(dfieldx_relative,this_linecut(:,1),this_linecut(:,2));
    dispyq = interp2(dfieldy_relative,this_linecut(:,1),this_linecut(:,2));
    this_line_disps = [dispxq,dispyq];
    linecutdisps(:,:,i) = this_line_disps;
end
mean_linecut_disps = mean(linecutdisps,3);
figure;
plot(mean_linecut_disps(:,1),mean_linecut_disps(:,2),'b-o');
hold on
plotFullDisplacementHexagons(gca);
% perp_unit_vec = [cos(perp_angle),sin(perp_angle)];

%% Set up Sine_Gordon
syms x w
delu = 2/pi*atan(exp(pi*x/w));
deluf = matlabFunction(delu);
xs = -5:0.01:5;
w = 4*ones(size(xs));
deluvals = deluf(w,xs);
figure
plot(xs,deluvals,'ro-');

% need to get projected data
[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
scale = hexagon_lattice_constant/sqrt(3);
zeropoint = v1-v2;
zeroed_disps = mean_linecut_disps - repmat(zeropoint,NinCut,1);
projected_disps = v2./(norm(v2))*zeroed_disps';
scaled_disps = projected_disps/scale;
% ascertain the distance
stepsize = 0.5;
xvals = linspace(0,2*r*stepsize,NinCut);
fdata = figure;
plot(fliplr(xvals),scaled_disps,'bo');
xlabel('Distance perpendicular to soliton boundary');
ylabel('Nondimensionalized stacking change');

%% Conduct fit
fitpred = @(w,translation,xs) deluf(w*ones(1,numel(xs)),xs-translation);
fitresid = @(params) fliplr(scaled_disps) - fitpred(params(1),params(2),xvals);
fitfun = @(params) rms(fitresid(params));
beta0 = [7,6];
options.Display = 'iter';
options.MaxFunEvals = 1e4;
beta = fminsearch(fitfun,beta0,options);
% Examine results
figure(fdata)
predvals = fitpred(beta(1),beta(2),xvals);
hold on
plot(xvals,predvals,'r-');
beta

%% Use the strain derivative function to get max strain across soliton boundary
units_strain_percent = @(w_soliton) 0.07104/w_soliton * 100;
strain_percent = units_strain_percent(beta(1))

