% ---------------------- Jacobian
% Mod NPK
function J = getjacobian(beta,fdiffstep,fun)

p = length(beta);

delta = zeros(size(beta));

yorg = fun(beta);

for k = 1:p

    if (beta(k) == 0)

        nb = sqrt(norm(beta));

        delta(k) = fdiffstep * (nb + (nb==0));

    else

        delta(k) = fdiffstep*beta(k);

    end

    yplus = fun(beta+delta);

    dy = yplus(:) - yorg(:);

%     dy(nans) = [];

    J(:,k) = dy/delta(k);

    delta(k) = 0;

end