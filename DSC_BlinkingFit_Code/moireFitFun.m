% fitting function, with ux, uy linear
function residuals = moireFitFun(ux,uy,x0s,y0s,thetam)
thetam = deg2rad(thetam);
b1 = ux - x0s*cos(thetam-1) + y0s*sin(thetam);
b2 = -uy + x0s*sin(thetam) + y0s*cos(thetam+1);
b = vertcat(b1,b2);
A1 = [1-cos(thetam),sin(thetam)];
A2 = [sin(thetam),cos(thetam)+1];
A = vertcat(repmat(A1,numel(b1),1),repmat(A2,numel(b2),1));
origin = A\b;  % completes separable linear least squares, now find residuals
xor = origin(1);
yor = origin(2);
uxpred = (x0s-xor)*(cos(thetam)-1) - (y0s-yor)*sin(thetam);
uypred = (x0s-xor)*sin(thetam) + (y0s-yor)*(cos(thetam)+1);
residuals = vertcat(ux,uy) - vertcat(uxpred,uypred);
end

