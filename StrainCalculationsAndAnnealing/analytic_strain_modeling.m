% analytic_strain_modeling.m

syms tr tm x y
rotmat_m1 = [cos(tm/2), -sin(tm/2); sin(tm/2), cos(tm/2)];
rotmat_m2 = [cos(tm/2), sin(tm/2); -sin(tm/2), cos(tm/2)];
trxy = exp(-x^2-y^2)*tr;
rotmat_r1 = [cos(trxy/2), -sin(trxy/2); sin(trxy/2), cos(trxy/2)];
rotmat_r2 = [cos(trxy/2), sin(trxy/2); -sin(trxy/2), cos(trxy/2)];

a1_0 = [x;y];
a2_0 = [x;y];

a1_m = rotmat_m1*a1_0;
a2_m = rotmat_m2*a2_0;
a1_r = rotmat_r1*a1_m;
a2_r = rotmat_r2*a2_m;

u_m = a1_m - a2_m;
disp('u_m(x,y)');
pretty(u_m);

exx_m = diff(u_m(1),'x');
exy_m = diff(u_m(1),'y');
eyx_m = diff(u_m(2),'x');
eyy_m = diff(u_m(2),'y');

u_r = a1_r - a2_r;
disp('u_r(x,y)');
pretty(u_r);

exx_obs_r = simplify(diff(u_r(1),'x'));
exy_obs_r = simplify(diff(u_r(1),'y'));
eyx_obs_r = simplify(diff(u_r(2),'x'));
eyy_obs_r = simplify(diff(u_r(2),'y'));

xbase = -2:0.01:2;
ybase = -2:0.01:2;
[xspace,yspace] = meshgrid(xbase,ybase);
exx_obs_fun = matlabFunction(exx_obs_r);
exy_obs_fun = matlabFunction(exy_obs_r);
eyx_obs_fun = matlabFunction(eyx_obs_r);
eyy_obs_fun = matlabFunction(eyy_obs_r);

theta_m = deg2rad(1.2);
theta_r = deg2rad(0.6);
figure;
imagesc(exx_obs_fun(theta_m,theta_r,xspace,yspace)); title('exx observed strain'); axis equal; set(gca,'yDir','normal'); colorbar;
figure;
imagesc(exy_obs_fun(theta_m,theta_r,xspace,yspace)); title('exy observed strain'); axis equal; set(gca,'yDir','normal'); colorbar;
figure;
imagesc(eyx_obs_fun(theta_m,theta_r,xspace,yspace)); title('eyx observed strain'); axis equal; set(gca,'yDir','normal'); colorbar;
figure;
imagesc(eyy_obs_fun(theta_m,theta_r,xspace,yspace)); title('eyy observed strain'); axis equal; set(gca,'yDir','normal'); colorbar;

exx_actual_1r = simplify(diff(a1_r(1),'x') - 1);  % would expect x and y components to vary linearly with themselves and themselves only. 
exy_actual_1r = simplify(diff(a1_r(1),'y'));
eyx_actual_1r = simplify(diff(a1_r(2),'x'));
eyy_actual_1r = simplify(diff(a1_r(2),'y') - 1);
exx_actual_fun = matlabFunction(exx_actual_1r);
exy_actual_fun = matlabFunction(exy_actual_1r);
eyx_actual_fun = matlabFunction(eyx_actual_1r);
eyy_actual_fun = matlabFunction(eyy_actual_1r);

figure;
imagesc(exx_actual_fun(theta_m,theta_r,xspace,yspace)); title('exx actual strain'); axis equal; set(gca,'yDir','normal'); colorbar;
figure;
imagesc(exy_actual_fun(theta_m,theta_r,xspace,yspace)); title('exy actual strain'); axis equal; set(gca,'yDir','normal'); colorbar;
figure;
imagesc(eyx_actual_fun(theta_m,theta_r,xspace,yspace)); title('eyx actual strain'); axis equal; set(gca,'yDir','normal'); colorbar;
figure;
imagesc(eyy_actual_fun(theta_m,theta_r,xspace,yspace)); title('eyy actual strain'); axis equal; set(gca,'yDir','normal'); colorbar;

% Would appear that the simplification of the analytic expression is beyond
% the simplification routine here. Try it numerically.
exx_ratio_ana = simplify(exx_obs_r/exx_actual_1r);
exy_ratio_ana = simplify(exy_obs_r/exy_actual_1r);
eyx_ratio_ana = simplify(eyx_obs_r/eyx_actual_1r);
eyy_ratio_ana = simplify(eyy_obs_r/eyy_actual_1r);

% No, it actually appears that the ratio is not exact, though close,
% especially at large strain induced by the twist.
exx_ratio = exx_obs_fun(theta_m,theta_r,xspace,yspace)./exx_actual_fun(theta_m,theta_r,xspace,yspace);
exx_ratiofun = matlabFunction(exx_ratio_ana);
figure; imagesc(exx_ratiofun(theta_m,theta_r,xspace,yspace)); caxis([1,3]); title('exx strain observed / exx in-plane strain'); axis equal; set(gca,'yDir','normal'); colorbar;
% If however the moire angle is ignored, then the observed strain appears
% to be a multiple of 2 of the actual interlayer strain. Note the reason
% this is interesting for exx is that it has no strain from the raw Moire,
% but the reconstruction process induces strain. 
figure; imagesc(exx_ratiofun(theta_m,theta_r,xspace,yspace)); caxis([1,3]); title('exx strain observed / exx in-plane strain'); axis equal; set(gca,'yDir','normal'); colorbar;
% Can we show this analytically?
exx_obs_nom = subs(exx_obs_r,tm,0);
exx_actual_nom = subs(exx_actual_1r,tm,0);
exx_nom_ratio = simplify(exx_obs_nom/exx_actual_nom);
pretty(exx_ratio_ana)
pretty(exx_nom_ratio)
pretty(exx_obs_r);
pretty(exx_actual_1r);
