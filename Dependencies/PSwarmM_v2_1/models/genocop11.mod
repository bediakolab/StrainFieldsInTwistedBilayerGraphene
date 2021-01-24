# Teste case 11 from Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures= Evolution Programs'
#   third edition 1996,Appendix C case 11 pp.223-240.

# aivaz 23-02-2008

var x{1..6}>=0;

minimize fx:
        6.5*x[1]-0.5*x[1]^2-x[2]-2*x[3]-3*x[4]-2*x[5]-x[6];

subject to c1:
    x[1]+2*x[2]+8*x[3]+x[4]+3*x[5]+5*x[6]-16<=0;

subject to c2:
    -8*x[1]-4*x[2]-2*x[3]+2*x[4]+4*x[5]-x[6]+1<=0;
    
subject to c3:
    2*x[1]+0.5*x[2]+0.2*x[3]-3*x[4]-x[5]-4*x[6]-24<=0;
    
subject to c4:
    0.2*x[1]+2*x[2]+0.1*x[3]-4*x[4]+2*x[5]+2*x[6]-12<=0;
  
subject to c5:
    -0.1*x[1]-0.5*x[2]+2*x[3]+5*x[4]-5*x[5]+3*x[6]-2<=0;  

subject to bound1 {j in 1..3}:
    x[j]<=10;
    
subject to bound2 {j in 4..5}:
    x[j]<=1;
    
subject to bound3:
    x[6]<=2;

# this problem have no initial guess (starting point) 
option reset_initial_guesses 1;

# global minimizer at [0 6 0 1 1 0] with fx = -11
