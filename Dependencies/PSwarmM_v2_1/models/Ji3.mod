# Problem 3 from Applied Mathematics and Computation
# Volume 185, Issue 1, 1 February 2007, Pages 382-387

# aivaz 23-02-2008

var x{1..2}>=0;

minimize fx:
    -((x[1]+3*x[2]+2)/(4*x[1]+x[2]+3)+
    (4*x[1]+3*x[2]+1)/(x[1]+x[2]+4));

subject to c1:
    -(x[1]+x[2])<=-1;


