function [ FX,FY ] = computeGradient( displacement_field, direction_string, permutation_number, xbound, ybound )
% direction_string should be either 'x' or 'y', which will give the
% component of the displacement field that we are differentiating.

TEST = 0;

uobj = DisplacementEquivalenceClass.empty;
utestobj = DisplacementEquivalenceClass.empty;
for q = 1:6400
    [a,b] = ind2sub([80,80],q);
    myobj = DisplacementEquivalenceClass(permute(displacement_field(a,b,:),[2,3,1]),direction_string,permutation_number,xbound,ybound);
    if TEST
        mytestobj = DisplacementEquivalenceClass(permute(displacement_field(a,b,:),[2,3,1]),'all',permutation_number,xbound,ybound);
        utestobj(q) = mytestobj;
    end
    uobj(q) = myobj;
end
uobj = reshape(uobj,[80,80]);
if TEST
    utestobj = reshape(utestobj,[80,80]);
end


h = 1:80;
n = 80;
f = uobj;
g = zeros(size(uobj));
for i = 1:80
    g(i,1) = (f(i,2) - f(i,1))/(h(2)-h(1));
    g(i,n) = (f(i,n) - f(i,n-1))/(h(end)-h(end-1));
    
    
end
% Take centered differences on interior points
h = h(3:n) - h(1:n-2);
for i = 1:80
    for j = 2:79
        g(i,j) = (f(i,j+1) - f(i,j-1)) ./ h(j-1);
        
        if i == 65 && j == 28
            if TEST
                utestobj(i,j+1).showDisplacement();
                utestobj(i,j-1).showDisplacement();
            end
        end
    end
    
end
FX = g;


% Do this again for the x direction
h = 1:80;
n = 80;
f = uobj';
g = zeros(size(uobj));
for i = 1:80
    g(i,1) = (f(i,2) - f(i,1))/(h(2)-h(1));
    g(i,n) = (f(i,n) - f(i,n-1))/(h(end)-h(end-1));
    
end
% Take centered differences on interior points
h = h(3:n) - h(1:n-2);
for i = 1:80
    for j = 2:79
        g(i,j) = (f(i,j+1) - f(i,j-1)) ./ h(j-1);
    end
end
FY = g';


end

