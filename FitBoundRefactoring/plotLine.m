function plotLine( figh, point1, point2, color, thickness )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n = 100;
figure(figh);
hold on;
xvals = linspace(point1(1),point2(1),n);
yvals = linspace(point1(2),point2(2),n);
if nargin < 4
    plot(xvals,yvals,'b-');
elseif nargin < 5
    plot(xvals,yvals,sprintf('%s-',color));
else
    plot(xvals,yvals,sprintf('%s-',color),'LineWidth',thickness);
end

end

