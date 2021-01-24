# Problem obtained from: Pintér, J.: Global optimization: software, test problems,
#    and applications. In: Handbook of Global Optimization, Vol. 2. Nonconvex
#    Optim. Appl., vol. 62. pp. 515–569. Kluwer Acad. Publ. Dordrecht (2002)
#
# Also available at: P. Parpas, B. Rustem and E.N. Pistikopoulos, Linearly constrained
#    global optimization and stochastic differential equations, J Glob Optim (2006)
#    DOI 10.1007/s10898-006-9026-z
#

option randseed 0;                             # to randomize

param n;                                       # externally set, number of variables
param m;                                       # externally set, number of linear constraints
param m_actives;                               # externally set, number of active linear constraints

param x_star{1..n} := Uniform(-10,10);         # randomly generated solution
param b{1..m};                                 # linear constraints independent term
param A{1..m,1..n} := Uniform(-10,10);         # linear constraints coefficient terms

param s;                                       # parameter
param g1;                                      # parameter
param g2;                                      # parameter

var x{1..n} <=10, >=-10;                       # variables

var P1 = sum{i in 1..n}((x[i]-x_star[i])^2)+sum{i in 1..n}((x[i]-x_star[i])^2);
var P2 = sum{i in 1..n}(x[i]-x_star[i]);

minimize fx:
    s*n*sum{i in 1..n}((x[i]-x_star[i])^2)+(sin(g1*P1))^2+(sin(g2*P2))^2;

subject to con {i in 1..m}:
    sum{j in 1..n} (x[j]*A[i,j])<=b[i];

data Pinter_n_m.dat;

data;

param s  := 0.025;
param g1 := 1;
param g2 := 1;

for {i in 1..m_actives} {
       let b[i] := sum{j in 1..n} (x_star[j]*A[i,j])+1e-3;                   # first m_actives constraints will be active
    }

for {i in (m_actives+1)..m} {
       let b[i] := sum{j in 1..n} (x_star[j]*A[i,j])+Uniform(1,10);     # last m-m_actives constraints will be inactive
    }

#display x_star;
#display A;
#display b;  
