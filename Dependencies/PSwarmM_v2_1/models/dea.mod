param d := 5;

param n := 100;

param a {1..d, 1..n}, default Uniform01();

var y;
var x {1..d} >= 0;

minimize ydist: (y - sum{j in 1..d} a[j,1]*x[j])^2 / sum{j in 1..d} x[j]*x[j];

subject to one_diside {i in 1..n}: sum {j in 1..d} a[j,i]*x[j] <= y;

