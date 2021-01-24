param N := 5;

var x{1..N};
param A{1..N};
param B{1..N};
param C{1..N};
param D{1..N};

minimize f:
-1*sum {i in 1..5} (A[i]*x[i] - 24345);

subject to cons1:
2.4*x[1] - x[2] >= 0;

subject to cons2:
-1.2*x[1] + x[2] >= 0;

subject to cons3:
60*x[1] - x[3] >= 0;

subject to cons4:
-20*x[1]+x[3] >= 0;

subject to cons5:
9.3*x[1]-x[4] >= 0;

subject to cons6:
-9.0*x[1]+x[4] >= 0;

subject to cons7:
7.0*x[1] - x[5] >= 0;

subject to cons8:
-6.5*x[1]+x[5] >= 0;

subject to cons9:
sum {i in 1..5} B[i]*x[i] >= 0;

subject to cons10:
sum {i in 1..5} C[i]*x[i] >= 0;

subject to cons11:
sum {i in 1..5} D[i]*x[i] >= 0;

subject to cons12:
-1*sum {i in 1..5} B[i]*x[i] + 294000 >= 0;

subject to cons13:
-1*sum {i in 1..5} C[i]*x[i] + 294000 >= 0;

subject to cons14:
-1*sum {i in 1..5} D[i]*x[i] + 294000 >= 0;

data;
var x:=
1   2.52
2   5.04
3   94.5
4   23.31
5   17.14;

param:
    A           B       C           D:=
1   -8720288.849        -145421.402 -155011.1084        -326669.5104
2   150512.5253     2931.1506   4360.53352      7390.68412
3   -156.6950325        -40.427932  12.9492344      -27.8986976
4   476470.3222     5106.192    10236.884       16643.076
5   729482.8271     15711.36    13176.786       30988.146;


