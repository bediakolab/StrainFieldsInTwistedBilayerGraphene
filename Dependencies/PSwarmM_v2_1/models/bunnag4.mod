#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 6
#

var x{1..6}, >=0;

minimize fx:
    -10.5*x[1]-7.5*x[2]-3.5*x[3]-2.547*x[4]-1.5*x[5]-10*x[6]-0.5*(sum{i in 1..5}(x[i]^2));

subject to con1:
    6*x[1]+3*x[3]+3*x[3]+2*x[4]+x[5]<=6.5;

subject to con2:
    10*x[1]+10*x[3]+x[6]<=20;

subject to bound1 {i in 1..5}: # x[6] not bounded above
    x[i]<=1;

# solution (0,1,0,1,1,20) f=-213
