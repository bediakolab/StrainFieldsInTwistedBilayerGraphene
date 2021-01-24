param N := 3;

var x{1..N};

minimize f:
-1*x[1]*x[2]*x[3];

subject to cons1:
1.8-x[1]-2*x[2] - 2*x[3] >= 0;

subject to cons2:
x[2] <= 1;

data;
var x:=
1   1
2   1
3   1;


