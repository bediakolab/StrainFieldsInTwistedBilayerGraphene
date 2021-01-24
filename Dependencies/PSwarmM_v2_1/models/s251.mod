param N := 3;

var x{1..N};

minimize f:
-1*x[1]*x[2]*x[3];

subject to cons1:
72 - x[1] - 2*x[2] - 2*x[3] >= 0;

subject to cons2{i in 1..N}:
0 <= x[i] <= 42;

data;
var x:=
1   10
2   10
3   10;

