function [f,grad_f] = QuadObjFun(x, H, g)

% It computes the function value and gradient 
% ---> for use in fmincom and DCA2

f = 0.5*x'*H*x+g'*x;
grad_f = H*x+g;
return;