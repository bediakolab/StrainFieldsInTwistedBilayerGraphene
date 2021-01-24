param N := 2;

var x{1..N};

minimize f:
100*(x[2]-x[1]^2)^2 + (1-x[1])^2;

subject to cons1:
x[1]/3 + x[2] + 0.1 >= 0;

subject to cons2:
-1*x[1]/3 + x[2] + 0.1 >= 0;

data;
var x:=
1   -1.2
2   1;


