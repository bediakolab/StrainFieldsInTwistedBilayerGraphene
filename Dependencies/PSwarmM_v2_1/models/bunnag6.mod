#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 8
#

param c{1..7};
var x{1..10}, >=0, <=1;

minimize fx:
    sum{i in 1..7} (c[i]*x[i])+10*sum{i in 8..10} (x[i])-5*sum{i in 1..7}(x[i]^2);

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
    
subject to con6:
    -7*x[1]-5*x[2]-2*x[3]-6*x[5]-6*x[6]-7*x[7]-6*x[8]+7*x[9]+7*x[10]<=-3;
    
subject to con7:
    x[1]-3*x[2]-3*x[3]-4*x[4]-x[5]-4*x[7]+6*x[9]<=1;
    
subject to con8:
    x[1]-2*x[2]+6*x[3]+9*x[4]-7*x[6]+9*x[7]-9*x[8]-6*x[9]+4*x[10]<=12;

subject to con9:
    -4*x[1]+6*x[2]+7*x[3]+2*x[4]+2*x[5]+6*x[7]+6*x[8]-7*x[9]+4*x[10]<=15;
    
subject to con10:
    sum{i in 1..10}(x[i])<=9;

subject to con11:
    -sum{i in 1..10}(x[i])<=-1;
    
data;

param c :=
1 -20
2 -80
3 -20
4 -50
5 -60
6 -90
7 0;


# solution (1,0.907, 0, 1, 0.715, 1, 0, 0.917, 1, 1) f=771.985
