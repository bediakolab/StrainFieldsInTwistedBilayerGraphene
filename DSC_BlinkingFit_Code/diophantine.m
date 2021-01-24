% diophantine.m
ourtheta = 1.2 * pi / 180;
% objfcn = @(m,n,theta) abs(acos(theta) - (m^2+4*m*n+n^2)/(2*(m^2+m*n+n^2)));
% to_use = @(invals) objfcn(invals(1),invals(2),ourtheta);
% fminsearch(to_use,[57,59])

objfcn = @(m,n) abs(ourtheta - acos((m^2+4*m*n+n^2)/(2*(m^2+m*n+n^2))));
fminsearch(to_use,[57,59]);
