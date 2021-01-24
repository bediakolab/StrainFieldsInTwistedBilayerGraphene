param N := 2;

var x{1..N};

minimize f:
(1/x[1])*log(2*log(x[2])/log(x[1]+x[2]));

subject to cons1:
1 - x[1] - x[2] >= 0;

subject to cons2{i in 1..2}:
x[i] >= 0.0001;

subject to cons3:
x[2] <= 1;

data;
var x:=
1   0.5
2   0.1;

