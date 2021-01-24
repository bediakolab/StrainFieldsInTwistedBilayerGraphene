# Teste case 7 from Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures= Evolution Programs'
#   third edition 1996,Appendix C case 11 pp.223-240.

# aivaz 23-02-2008

var x{1..6}>=0;

minimize fx:
        -10.5*x[1]-7.5*x[2]-3.5*x[3]-2.5*x[4]-1.5*x[5]-10*x[6]-0.5*sum{i in 1..5}(x[i]^2);

subject to c1:
    6*x[1]+3*x[2]+3*x[3]+2*x[4]+x[5]-6.5<=0;

subject to c2:
    10*x[1]+10*x[3]+x[6]-20<=0;

subject to bound1 {j in 1..5}:
    x[j]<=1;

subject to bound2:
    x[6]<=100;

# this problem have no initial guess (starting point) 
option reset_initial_guesses 1;

# global minimizer at [0 1 0 1  1 20] with fx = -213
