param N := 4;

var x{1..N};

minimize f:
(x[1]+10*x[2])^2+5*(x[3]-x[4])^2+(x[2]-2*x[3])^4+10*(x[1]-x[4])^4;

subject to cons1:
x[1]+x[2]+x[3]+x[4] - 1 >= 0;

subject to cons2{i in 1..N}:
x[i] <= 20;
data;
var x:=
1   3
2   -1
3   0
4   1;


