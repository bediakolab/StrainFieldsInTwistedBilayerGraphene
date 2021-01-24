function [ pred_vals, J ] = HamishTrigFittingFunction( displacement_guesses, prefactors )
% Objective function for interferometry fitting of the 4DSTEM data. Uses
% Hamish's approach of not only fitting the displacement vectors, but also
% the magnitudes of the cosine prefactors. 
%
% Because this function needs to be quite fast, there will be no nan
% handling of any sort. Optimizations are allowed to run outside of the
% half-hexagon displacement bounds, and we will worry about the
% re-conversion after the fact.
%
% Inputs: displacement_guesses is an 2xn array. x and y components will
% have been flattened out and concatenated horizontally in caller before
% passing to this fitting function. Prefactors is a 12x1 vector.
%
% Multiple inputs are ok here because the optimizer never sees this
% function directly.
%
% Nathanael Kazmierczak, 04/06/2020

mat = [  % Lattice plane normal vectors
   0.866025403784439  -0.500000000000000
   0.000000000000000  -1.000000000000000
  -0.866025403784439  -0.500000000000000
   0.866025403784439  -0.500000000000000
   0.000000000000000  -1.000000000000000
  -0.866025403784438  -0.500000000000000
   1.000000000000000   0.000000000000000
   0.500000000000000  -0.866025403784439
  -0.500000000000000  -0.866025403784439
   1.000000000000000   0.000000000000000
   0.500000000000000  -0.866025403784438
  -0.500000000000000  -0.866025403784439];

plane_sep = [ 
   2.131288518713503
   2.131288518713503
   2.131288518713503
   2.131288518713503
   2.131288518713503
   2.131288518713503
   1.230500000000000
   1.230500000000000
   1.230500000000000
   1.230500000000000
   1.230500000000000
   1.230500000000000];

mags = mat*displacement_guesses;  % produces nx12 matrix of projection magntudes for 12 disks
mul = pi./plane_sep;
nondim_predvals = cos(mags.*mul).^2;
pred_vals = prefactors.*nondim_predvals;
pred_vals = pred_vals(:);  % stacks 12-vectors on top of each other vertically.

%% Jacobian calculation
% Added 03/27/2020 after meeting with Colin and Hamish at NCEM.
% For this basic case, it is a 12x2 Jacobian.
% NPK 04/06/2020: letting the image pixel numel be k, the new full dimension is
% 12k by 12 + 2*k. But that's waaay too big to construct explicitly, so
% we will do it sparsely.
if nargout > 1
    k = size(displacement_guesses,2);  % number of image pixels
    % Assemble prefactor portion of the Jacobian
    i_prefactor = (1:(k*12))';
    j_prefactor = mod(i_prefactor-1,12)+1;
    v_prefactor = -nondim_predvals(:);
    
    % Assemble displacement portion of the Jacobian
    mm = mags.*mul;
    base_portion = prefactors.*sin(2*mm);
    xpor = mul.*mat(:,1).*base_portion;
    ypor = mul.*mat(:,2).*base_portion;
    bo = ones(12,k);
    xpor_jidx = bo.*((1:2:(2*k-1)) + 12);
    ypor_jidx = bo.*((2:2:(2*k)) + 12);
    
%     Jold = zeros(12,2);
%      Jold(:,1) = 2*prefactors.*cos(mul.*mags(:,2)).*sin(mul.*mags(:,2)).*mul.*mat(:,1);
%      Jold(:,2) = 2*prefactors.*cos(mul.*mags(:,2)).*sin(mul.*mags(:,2)).*mul.*mat(:,2);

    % Final call to construct Jacobian
    i = vertcat(i_prefactor,i_prefactor,i_prefactor);  % because it's the same linear i index for each one.
    j = vertcat(j_prefactor,xpor_jidx(:),ypor_jidx(:));
    v = vertcat(v_prefactor,xpor(:),ypor(:));
    J = sparse(i,j,v);
end

end

