# Teste case #3 from Michalewicz, Z., Evolutionary Computation Techniques
#   for Nonlinear Programming Problems, International Transactions in
#   Operational Research, Vol.1, No.2, 1994, pp.223-240.

# Teste case 11 from Michalewicz ,Zbigniew 'Genetic Algorithms+ Data Structures= Evolution Programs'
#   third edition 1996,Appendix C case 11 pp.223-240.

# aivaz 23-02-2008

var x{1..2}>=0;

minimize fx:
    (if (x[1]<2) then (x[2]+1e-5*(x[2]-x[1])^2-1.0) else 
        (if (x[1]<4) then (1/(27*sqrt(3))*((x[1]-3)^2-9)*x[2]^3) else
            (1/3*(x[1]-2)^3+x[2]-11/3)));

subject to c1:
    x[1]/sqrt(3)-x[2]>=0;

subject to c2:
    -x[1]-sqrt(3)*x[2]+6>=0;

subject to bound:
    x[1]<=6;
    
# this problem have no initial guess (starting point)
option reset_initial_guesses 1;

# three global minimizers at (0,0) (3,sqrt(3)) and (4,0) with fx = -1
