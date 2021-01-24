#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 1
#

var x{1..3}, >=0, <=3;

minimize fx:
    9-8*x[1]-6*x[2]-4*x[3]+2*x[1]^2+2*x[2]^2+x[3]^2+2*x[1]*x[2]+2*x[1]*x[3];

subject to con:
    x[1]+x[2]+2*x[3]<=3;
