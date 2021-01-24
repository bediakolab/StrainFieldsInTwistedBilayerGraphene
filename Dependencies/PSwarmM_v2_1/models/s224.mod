param N := 2;

var x{1..N};

minimize f:
2*x[1]^2 + x[2]^2 - 48*x[1] - 40*x[2];

subject to cons1:
x[1] + 3*x[2] >= 0;

subject to cons2:
18 - x[1] - 3*x[2] >= 0;

subject to cons3:
x[1] + x[2] >= 0;

subject to cons4:
8 - x[1] - x[2] >= 0;

subject to cons5:
0 <= x[1] <= 6;

subject to cons6:
0 <= x[2] <= 6;

data;
var x:=
1   0.1
2   0.1;


