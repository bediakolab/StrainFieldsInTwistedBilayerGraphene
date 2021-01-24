# This version works

#      1                 2
#     /       / 1 + y'(x)  \
# min |  sqrt | ---------- |  dx
#     /       \    y(x)    /
#     0
#

#param n := 512;
param n := 128; # aivaz 26-02-2008

param x {j in 0..n} := (j/n)^4;

var y {j in 0..n} >= 0;
var dydx {j in 1..n} = (y[j] - y[j-1])/(x[j] - x[j-1]);
#var f {j in 1..n} = sqrt( (1+dydx[j]^2)/((y[j]+y[j-1])/2) );
var f {j in 1..n} = sqrt( (1+dydx[j]^2)/y[j-1] );

minimize time: sum {j in 1..n} f[j]*(x[j]-x[j-1]) ;

subject to y0: y[0] = 1.0e-12;
subject to yn: y[n] = 1;

# doesn't work without monotonicity assumption
subject to monotone {j in 1..n}: dydx[j] >= 0;

let {j in 0..n} y[j] := x[j]^0.7;



