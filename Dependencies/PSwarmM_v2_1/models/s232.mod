param N := 2;

var x{1..N};

minimize f:
-1*(9-(x[1] - 3)^2)*(x[2]^3/(27*sqrt(3)));

subject to cons1:
x[1]/sqrt(3) - x[2] >= 0;

subject to cons2:
x[1] + sqrt(3)*x[2] >= 0;

subject to cons3:
6 - x[1] - sqrt(3)*x[2] >= 0;

subject to cons4:
x[1] >= 0;

subject to cons5:
x[2] >= 0;

data;
var x:=
1   2
2   0.5;


