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

#   Source:  problem 2 in A. Zecevic, "Contribution to methods 
#            of external penalty functions - algorithm MNBS"
#            Advanced Business School, Belgrade, 
#            (whatever is left of) Yugoslavia.

#   SIF input: Nick Gould, April 1993.

#   classification QLR2-AN-2-2

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
    2*x[2]^2-2*x[1]-3*x[2];
subject to cons1:
    x[1]+x[2] <= 2.0;
subject to cons2:
    x[1]+4*x[2] <= 4.0;
subject to cons3:
    0 <= x[1] <= 10;
subject to cons4:
    0 <= x[2] <= 10;

data;
param xinit:= 1 0.1 2 -0.1;

