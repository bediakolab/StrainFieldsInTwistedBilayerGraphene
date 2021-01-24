param N := 3;

var x{1..N} >= 0;
param a{1..8,1..3};
param c{1..8}:= 1;

minimize f:
sum {j in 1..8} c[j]*(sum {i in 1..3} (a[j,i] - x[i])^2)^0.5;

subject to cons1:
30 - 3*x[1] - 3*x[3] >= 0;

data;
var x:=
1   0
2   2
3   0;

param a:
    1   2   3:=
1   0   0   0
2   10  0   0
3   10  10  0
4   0   10  0
5   0   0   10
6   10  0   10
7   10  10  10
8   0   10  10;

    

