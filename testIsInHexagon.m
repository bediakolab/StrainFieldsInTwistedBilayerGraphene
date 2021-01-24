% testIsInHexagon.m

x0 = 0; y0 = 0; d = 1;
latticepoints = getHexagon(x0,y0,d);
[xsort] = sortrows(latticepoints);
finalsort = [xsort(1:2,:); xsort(4,:); xsort(6,:); xsort(5,:); xsort(3,:); xsort(1,:)];
figure;
plot(finalsort(:,1),finalsort(:,2));
while 1
    [x,y] = ginput(1);
    tf = isInHexagon(x,y,x0,y0,d)
end

