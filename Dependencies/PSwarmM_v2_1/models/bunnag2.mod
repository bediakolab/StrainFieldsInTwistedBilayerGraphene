#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 3
#

var x{1..4}, >=0, <=4;

minimize fx:
    x[1]^0.6+2*x[2]^0.6-2*x[2]+2*x[3]-x[4];

subject to con1:
    x[1]+2*x[3]<=4;

subject to con2:
    -3*x[1]+x[4]<=1;


# solution (4/3, 4, 0, 0) f=-2.07
