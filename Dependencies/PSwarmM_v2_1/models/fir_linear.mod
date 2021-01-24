param pi := 4*atan(1);

param beta := 0.01;

param omega_stop := 2*pi/3;
param omega_pass :=   pi/2;
#param step := pi/180;
param step := pi/45; # aivaz 26-02-2008
set OMEGA_STOP := {omega_stop..pi by step};
set OMEGA_PASS := {0..omega_pass  by step};

param n := 20;  # must be even
param n2 := n/2;

var h {0..n2-1};
var t;

minimize spread: t;

subject to passband_up_bnds {omega in OMEGA_PASS}:
           2* sum {k in 0..n2-1} h[k]*cos((k-(n-1)/2)*omega) <= 1+t;

subject to passband_lo_bnds {omega in OMEGA_PASS}:
    1-t <= 2* sum {k in 0..n2-1} h[k]*cos((k-(n-1)/2)*omega);

subject to stopband_bnds {omega in OMEGA_STOP}:
    -beta <= 2* sum {k in 0..n2-1} h[k]*cos((k-(n-1)/2)*omega) <= beta;

