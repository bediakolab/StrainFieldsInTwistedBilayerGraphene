#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 9
#

param c{1..10};
var x{1..10}, >=0, <=1;

minimize fx:
    sum{i in 1..10} (c[i]*x[i])-50*sum{i in 1..10}(x[i]^2);

subject to con1:
    -2*x[1]-6*x[2]-x[3]-3*x[5]-3*x[6]-2*x[7]-6*x[8]-2*x[9]-2*x[10]<=-4;

subject to con2:
    6*x[1]-5*x[2]+8*x[3]-3*x[4]+x[6]+3*x[7]+8*x[8]+9*x[9]-3*x[10]<=22;

subject to con3:
    -5*x[1]+6*x[2]+5*x[3]+3*x[4]+8*x[5]-8*x[6]+9*x[7]+2*x[8]-9*x[10]<=-6;

subject to con4:
    9*x[1]+5*x[2]-9*x[4]+x[5]-8*x[6]+3*x[7]-9*x[8]-9*x[9]-3*x[10]<=-23;

subject to con5:
    -8*x[1]+7*x[2]-4*x[3]-5*x[4]-9*x[5]+x[6]-7*x[7]-x[8]+3*x[9]-2*x[10]<=-12;
    
    
data;

param c :=
1 48
2 42
3 48
4 45
5 44
6 41
7 47
8 42
9 45
10 46;


# solution (1, 0, 0, 1, 1, 1, 0, 1, 1, 1) f=-39

