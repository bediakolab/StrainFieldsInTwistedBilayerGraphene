classdef vectorLine
    % vectorLine.m
    % Because NPK doesn't like doing this by hand each time.
    %
    % If slope is vertical, then will use x predictions instead.
    
    properties
        vec
        point
        m
    end
    
    methods
        function obj = vectorLine(vec,point)
            obj.vec = vec;
            obj.point = point;
            obj.m = vec(2)/vec(1);
        end
        
        function tf = lt(obj,otherpoint)
            y = otherpoint(:,2);
            x = otherpoint(:,1);
            if abs(obj.m) < 1000
                ypred = obj.getypred(otherpoint);
                tf = ypred < y;
            else
                tf = obj.point(1) < x;
            end
        end
        
        function tf = gt(obj,otherpoint)
            y = otherpoint(:,2);
            x = otherpoint(:,1);
            if abs(obj.m) < 1000
                ypred = obj.getypred(otherpoint);
                tf = ypred > y;
            else
                tf = obj.point(1) > x;
            end
        end
        
        function tf = le(obj,otherpoint)
            y = otherpoint(:,2);
            x = otherpoint(:,1);
            if abs(obj.m) < 1000
                ypred = obj.getypred(otherpoint);
                tf = ypred <= y;
            else
                tf = obj.point(1) <= x;
            end
        end
        
        function tf = ge(obj,otherpoint)
            y = otherpoint(:,2);
            x = otherpoint(:,1);
            if abs(obj.m) < 1000
                ypred = obj.getypred(otherpoint);
                tf = ypred >= y;
            else
                tf = obj.point(1) >= x;
            end
        end
        
        function ypred = getypred(obj,otherpoint)
            x = otherpoint(:,1);
            y = otherpoint(:,2);
            x0 = obj.point(1);
            y0 = obj.point(2);
            ypred = obj.m*(x - x0) + y0;
        end
        
        
    end
    
end

