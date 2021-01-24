# Problem 1 from Applied Mathematics and Computation 
# Volume 185, Issue 1, 1 February 2007, Pages 382-387 

# aivaz 23-02-2008

var x{1..3}>=0;

minimize fx:
    -((4*x[1]+3*x[2]+3*x[3]+50)/(3*x[2]+3*x[3]+50)+
    (3*x[1]+4*x[3]+50)/(4*x[1]+4*x[2]+5*x[3]+50)+
    (x[1]+2*x[2]+5*x[3]+50)/(x[1]+5*x[2]+5*x[3]+50)+
    (x[1]+2*x[2]+4*x[3]+50)/(5*x[2]+4*x[3]+50));
    
subject to c1:
    2*x[1]+x[2]+5*x[3]<=10;
    
subject to c2:
    x[1]+6*x[2]+3*x[3]<=10;
    
subject to c3:
    5*x[1]+9*x[2]+2*x[3]<=10;
    
subject to c4:
    9*x[1]+7*x[2]+3*x[3]<=10;
