% fitting_schematic.m

% [ pred_vals ] = trigFittingFunctions( DSC_guess, scaling_constant, bound_handling_flag );
xbase = -1.3:0.01:1.3;
ybase = 0:0.01:1.5;
[xdisp,ydisp] = meshgrid(xbase,ybase);
dispfit = [xdisp(:),ydisp(:)]';
prefactors = ones(12,1);
[ pred_vals ] = HamishTrigFittingFunction( dispfit, prefactors );
k = numel(pred_vals);
n = numel(xdisp);
rp = reshape(pred_vals,[12,n]);


testdisps = dispfit';
[ rz_disps ] = extendedZoneDisp2ReducedZoneDisp( testdisps );
tol = 1e-10;
tf1 = abs(testdisps(:,1) - rz_disps(:,1)) < tol;
tf2 = abs(testdisps(:,2) - rz_disps(:,2)) < tol;
tf = tf1 & tf2;
tf = reshape(tf,size(xdisp));

for i = [1:3,7:9]
    di = reshape(rp(i,:)',size(xdisp));
    di(~tf) = nan;
    figure;
    h = contourf(xbase,ybase,di,1000,'LineStyle','none');
    shading flat;
    hold on; plotFullDisplacementHexagons(gca);
    title(sprintf('Predicted relative intensities for disk %d',i));
    colormap(fire(100));
%     colormap(gray(1000));
    colorbar
    xlim([-1.3,1.3]);
    ylim([-0.3,1.7]);
    caxis([0,1]);
    xlabel('x displacement (Angstrom)');
    ylabel('y displacement (Angstrom)');
end

% try cam angle

% % Trying Walter Robinson's stacking code:
% %https://www.mathworks.com/matlabcentral/answers/222993-how-to-stack-layers-of-2d-plot-in-horizontally-instead-of-vertically
% number_of_layers = 3;
% X_spacing = 1.234;     %adjust per tastes, 1 is valid
% first_X = 0.5;         %adjust per tastes, 1 is valid
% x = -2:.1:2;
% y = x;
% nx = length(x);
% ny = length(y);cam an
% z = rand(nx,ny,number_of_layers);
% [Ycoord, Zcoord] = ndgrid(x, y);
% for i = 1 : number_of_layers
%     di = reshape(rp(i,:)',size(xdisp));
%     di(~tf) = nan;
%     figure;
%     pcolor(xbase,ybase,di);
%     shading flat;
%     hold on; plotFullDisplacementHexagons(gca);
%     title(sprintf('Disk %d',i));
%     colormap(gray(10000));
%     xlim([-1.3,1.3]);
%     ylim([-0.3,1.7]);
%     
%     
%   Xcoord = (first_X + X_spacing * (i-1)) * ones(nx, ny);
%   surf(Xcoord, Ycoord, Zcoord, z(:,:,i));
%   hold on
% end
% ylim([-2.5 2.5]);
% zlim([-2.5 2.5]);

