# Teste case 9 from Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures= Evolution Programs'
#   third edition 1996,Appendix C case 11 pp.223-240.

# aivaz 23-02-2008

var x{1..3}>=0;

minimize fx:
        -((3*x[1]+x[2]-2*x[3]+0.8)/(2*x[1]-x[2]+x[3])+(4*x[1]-2*x[2]+x[3])/(7*x[1]+3*x[2]-x[3]));

subject to c1:
    x[1]+x[2]-x[3]-1<=0;

subject to c2:
    -x[1]+x[2]-x[3]+1<=0;
    
subject to c3:
    12*x[1]+5*x[2]+12*x[2]-34.8<=0;
    
subject to c4:
    12*x[1]+12*x[2]+7*x[3]-29.1<=0;
    
subject to c5:
    -6*x[1]+x[2]+x[3]+4.1<=0;

subject to bound {j in 1..3}:
    x[j]<=10;


# this problem have no initial guess (starting point) option
option reset_initial_guesses 1;

# global minimizer at [1 0 0] with fx = -2.4714
