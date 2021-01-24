function [ predvals,J ] = ellipticGaussianPredfun( space,x0,y0,sigmax,sigmay, A, B, C )
%Note that there is no covariance permitted here, but only a rotation
%parameter that rotates the coordinate basis.
%
% Nathanael Kazmierczak, 04/01/2020

% space, the coordinates, needs to be an nx2 array with x in first col and
% y in second
rotmat = [cos(C), -sin(C); sin(C), cos(C)];
us = repmat([x0,y0],size(space,1),1);
newspace = (rotmat*(space-us)')'+us;
x = newspace(:,1);
y = newspace(:,2);

% @(x0,y0,sigma1,sigma2,A,B,C) A - B*mvnpdf(space,[x0,y0],[]*[sigma1,0;0,sigma2]);
xp = (x-x0).^2./(2*sigmax.^2);
yp = (y-y0).^2./(2*sigmay.^2);
predvals = A - B*exp(-(xp + yp));

% Compute analytic Jacobian
if nargout > 1
    % order of parameters is [x0,y0,sigma1,sigma2, A, B, C]
    q = predvals - A;
    gx = x;
    gy = y;
    qsx = -(gx-x0)./(sigmax.^2);
    qsy = -(gy-y0)./(sigmay.^2);
    dgxt = -(x-x0).*sin(C)-(y-y0).*cos(C);
    dgyt = (x-x0).*cos(C)-(y-y0).*sin(C);
    
    dg2xd0 = -2*(x-x0).^2.*cos(C).*sin(C) - 2*(x-x0).*(y-y0).*(cos(C).^2-sin(C).^2) + 2*(y-y0).^2.*sin(C).*cos(C);
    dg2yd0 = -2*(y-y0).^2.*cos(C).*sin(C) + 2*(x-x0).*(y-y0).*(cos(C).^2-sin(C).^2) + 2*(x-x0).^2.*sin(C).*cos(C);
    
    
    
    J(:,1) = q.*(qsx.*(-cos(C)) + qsy.*(-sin(C)));
    J(:,2) = q.*(qsx.*(sin(C)) + qsy.*(-cos(C)));
    J(:,3) = q.*(gx-x0).^2./(sigmax.^3);
    J(:,4) = q.*(gy-y0).^2./(sigmay.^3);
    J(:,5) = 1;
    J(:,6) = q./B;
%     J(:,7) = q.*(qsx.*dgxt + qsy.*dgyt);
% x0,y0,sigmax,sigmay, A, B, C
    J(:,7) = q.*(-dg2xd0./(2*sigmax.^2) - dg2yd0./(2*sigmay.^2));
end



end

