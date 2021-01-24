# Problem g01 from: TP Runarsson, X Yao 'Stochastic Ranking for Constrained Evolutionary Optimization'
# IEEE TRANSACTIONS ON EVOLUTIONARY COMPUTATION, 2000 

# aivaz 23-02-2008

var x {1..13};

minimize fx: 5*sum{j in 1..4} x[j] - 5*sum{j in 1..4} x[j]^2 - sum{j in 5..13} x[j];

subject to c1:
    2*x[1]+2*x[2]+x[10]+x[11]-10<=0;
    
subject to c2:
    2*x[1]+2*x[3]+x[10]+x[12]-10<=0;
    
subject to c3:
    2*x[2]+2*x[3]+x[11]+x[12]-10<=0;
    
subject to c4:
    -8*x[1]+x[10]<=0;
    
subject to c5:
    -8*x[2]+x[11]<=0;
    
subject to c6:
    -8*x[3]+x[12]<=0;
    
subject to c7:
    -2*x[4]-x[5]+x[10]<=0;
    
subject to c8:
    -2*x[6]-x[7]+x[11]<=0;
    
subject to c9:
    -2*x[8]-x[9]+x[12]<=0;

subject to bounds1 {j in 1..9}:
    0<=x[j]<=1;
    
subject to bounds2 {j in 10..12}:
    0<=x[j]<=100;
    
subject to bounds3:
    0<=x[13]<=1;
