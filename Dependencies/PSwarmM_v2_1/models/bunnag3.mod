#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 5
#

var x{1..5}, >=0;

minimize fx:
    x[1]^0.6+x[2]^0.6+x[3]^0.6-4*x[3]-2*x[4]+5*x[5];

subject to con1:
    x[1]+2*x[4]<=4;

subject to con2:
    3*x[1]+3*x[4]+x[5]<=4;
    
subject to con3:
    2*x[2]+4*x[4]+2*x[5]<=6;

subject to bound1:
    x[1]<=3;
subject to bound2:
    x[3]<=4;
subject to bound3:
    x[5]<=2;
    
# solution (0.67,2,4,0,0) f=-11.96

