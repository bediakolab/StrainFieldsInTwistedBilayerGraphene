function [annulus_mask] = isInAnnulus(xspace,yspace,x0,y0,rinner,router)
% rinner should be greater than router
assert(rinner < router);
maskinner = isInCircle(xspace,yspace,x0,y0,rinner);
maskouter = isInCircle(xspace,yspace,x0,y0,router);
annulus_mask = maskouter & ~maskinner;

end

