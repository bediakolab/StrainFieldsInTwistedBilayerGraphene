function [f]=ampl_obj(x)
% Objective function should be prepared to receive a set of columnwise points

m=size(x,2);
f=zeros(m,1);
for i=1:m
    f(i) = matampl(x(:,i));
end

return;