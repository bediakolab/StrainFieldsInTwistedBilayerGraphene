function [f,grad_f] = RBFObjFun(x, Y, lambda, n, np)

% It computes the function value and gradient 
% ---> for use in fmincom and DCA2

f = 0;
grad_f = zeros(n,1);
for i=1:np
    v = x-Y(:,i);
    nv = norm(v);
    f = f + lambda(i)*nv^3;
    grad_f = grad_f + lambda(i)*3*nv*v;
end

f = f + lambda(np+1) + lambda(np+2:np+n+1)'*x;
grad_f = grad_f + lambda(np+2:np+n+1);
return;