# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source:
#   A.R. Conn, N. Gould and Ph.L. Toint,
#   "The LANCELOT User's Manual",
#   Dept of Maths, FUNDP, 1991.

#   SIF input: Ph. Toint, Jan 1991.

#   classification SLR2-AN-2-1

var a>=0;
var b;
param c := 0.85;
param x{1..5};
param y{1..5};

minimize f:
    0.5*sum {i in 1..5} (a*x[i]+b-y[i])^2;
subject to cons1:
    a+b <= c;

data;
param x:=
1   0.1
2   0.3
3   0.5
4   0.7
5   0.9;

param y:=
1   0.25
2   0.3
3   0.625
4   0.701
5   1.0;

