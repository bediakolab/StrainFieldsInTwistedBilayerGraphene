param N := 5;

var x{1..N} := 1;
param D{1..6,1..5};
param DtD{i in 1..5, j in 1..5} := sum {k in 1..6} D[k,i]*D[k,j];
param d{1..6};

minimize f:
sum {i in 1..5} x[i]*(sum {j in 1..5} DtD[i,j]*x[j]) - 2*sum {i in 1..6} d[i]*(sum {j in 1..5} D[i,j]*x[j]) + sum 
{i in 1..6} d[i]^2;

subject to cons1:
-1*x[1] - x[2] - x[3] - x[4] - x[5] + 5 >= 0;

subject to cons2:
10*x[1] + 10*x[2] - 3*x[3] + 5*x[4] + 4*x[5] - 20 >= 0;

subject to cons3:
-8*x[1] + x[2] - 2*x[3] - 5*x[4] + 3*x[5] + 40 >= 0;

subject to cons4:
8*x[1] - x[2] + 2*x[3] + 5*x[4] - 3*x[5] - 11 >= 0;

subject to cons5:
-4*x[1] - 2*x[2] + 3*x[3] - 5*x[4] + x[5] + 30 >= 0;

data;
param D:
    1   2   3   4   5:=
1   -74 80  18  -11 -4
2   14  -69 21  28  0
3   66  -72 -5  7   1
4   -12 66  -30 -23 3
5   3   8   -7  -4  1
6   4   -12 4   4   0;

param d:=
1   51
2   -61
3   -56
4   69
5   10
6   -12;

