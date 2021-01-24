# Objective: linear
# Constraints: linear

#param m default 14;    # must be even
#param n default 21;
param m default  6; # must be even
param n default 12;

set X := {0..n};
set Y := {0..m};

set NODES := X cross Y;     # A lattice of Nodes

set ANCHORS within NODES := { x in X, y in Y : 
                  x == 0 && y >= floor(m/3) && y <= m-floor(m/3) };

param load {(x,y) in NODES, d in 1..2} default 0;

param gcd {x in -n..n, y in -n..n} :=
    (if x <  0 then gcd[-x,y] else 
    (if x == 0 then y else 
    (if y < x  then gcd[y,x] else 
    (gcd[y mod x, x])
    )
    )
    );

set ARCS := { (xi,yi) in NODES, (xj,yj) in NODES: 
    abs( xj-xi ) <= 3 && abs(yj-yi) <=3 && abs(gcd[ xj-xi, yj-yi ]) == 1 
    && ( xi > xj || (xi == xj && yi > yj) )
    };

param length {(xi,yi,xj,yj) in ARCS} := sqrt( (xj-xi)^2 + (yj-yi)^2 );

var w {NODES, 1..2};

maximize work: 
    sum {(x,y) in NODES, d in 1..2} load[x,y,d]*w[x,y,d];
#minimize work: 
#    - sum {(x,y) in NODES, d in 1..2} load[x,y,d]*w[x,y,d];

subject to element_eqs_p {(xi,yi,xj,yj) in ARCS}: 
    (xj-xi)/length[xi,yi,xj,yj] * (w[xi,yi,1] - w[xj,yj,1])
    +
    (yj-yi)/length[xi,yi,xj,yj] * (w[xi,yi,2] - w[xj,yj,2])
    <= length[xi,yi,xj,yj];

subject to element_eqs_m {(xi,yi,xj,yj) in ARCS}: 
    -length[xi,yi,xj,yj]
    <= 
    (xj-xi)/length[xi,yi,xj,yj] * (w[xi,yi,1] - w[xj,yj,1])
    +
    (yj-yi)/length[xi,yi,xj,yj] * (w[xi,yi,2] - w[xj,yj,2])
    ;

subject to anchor {(x,y) in ANCHORS, d in 1..2}:
    w[x,y,d] = 0;

let load[n,m/2,2] := -1;

let {(x,y) in NODES, d in 1..2} w[x,y,d] := 0.001;

