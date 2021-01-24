#
# D. Bunnag and M. Sun, "Genetic algorithm for constrained global optimization in continuous variables",
# Applied Mathematics and Computation, 171 (2005) 604-636
#
# Problem 7
#

var x{1..6}, >=0;

minimize fx:
    6.5*x[1]-0.5*x[1]^2-x[2]-2*x[3]-3*x[4]-2*x[5]-x[6];

subject to con1:
    x[1]+2*x[2]+8*x[3]+x[4]+3*x[5]+5*x[6]<=16;

subject to con2:
    -8*x[1]-4*x[2]-2*x[3]+2*x[4]+4*x[5]-x[6]<=-1;

subject to con3:
    2*x[1]+0.5*x[2]+0.2*x[3]-3*x[4]-x[5]-4*x[6]<=24;

subject to con4:
    0.2*x[1]+2*x[2]+0.1*x[3]-4*x[4]+2*x[5]+2*x[6]<=12;

subject to con5:
    -0.1*x[1]-0.5*x[2]+2*x[3]+5*x[4]-5*x[5]+3*x[6]<=3;

subject to bound1:
    x[1]<=2;

subject to bound2:
    x[2]<=8;

subject to bound3:
    x[3]<=2;

subject to bound4:
    x[4]<=1;

subject to bound5:
    x[5]<=1;

subject to bound6:
    x[6]<=2;
    
# solution (0,6,0,1,1,0) f=-11.005
