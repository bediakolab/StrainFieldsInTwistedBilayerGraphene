# Teste case 10 from Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures= Evolution Programs'
#   third edition 1996,Appendix C case 11 pp.223-240.

# aivaz 23-02-2008

var x{1..4}>=0;

minimize fx:
        -(x[1]^0.6+x[2]^0.6-6*x[1]-4*x[3]+3*x[4]);

subject to c1:
    x[1]+2*x[3]-4<=0;

subject to c2:
    x[2]+2*x[4]-4<=0;

subject to bound1:
    x[1]<=3;

subject to bound2 {j in 2..3}:
    x[j]<=10;

subject to bound3:
    x[4]<=1;

# this problem have no initial guess (starting point)
option reset_initial_guesses 1;

# global minimizer at [4/3 4 0 0] with fx = -4.5142 ??
