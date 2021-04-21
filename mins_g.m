function g = mins_g(x)

n = size(x,1)-1;
h = 1/n;
a = x(1:end-1,1:end-1)-x(1:end-1,2:end);
b = x(2:end,2:end)-x(1:end-1,2:end);
c = x(2:end,2:end)-x(2:end,1:end-1); 
d = x(1:end-1,1:end-1)-x(2:end,1:end-1); 
C1 = sqrt(h^2 + a.^2 + b.^2); 
C2 = sqrt(h^2 + c.^2 + d.^2); 
g = zeros(n+1,n+1);
g(2:end-1,2:end-1) = a(2:end,2:end)./C1(2:end,2:end) + d(2:end,2:end)./C2(2:end,2:end) + ...
                     b(1:end-1,1:end-1)./C1(1:end-1,1:end-1) + c(1:end-1,1:end-1)./C2(1:end-1,1:end-1) - ...
                     (a(2:end,1:end-1)+b(2:end,1:end-1))./C1(2:end,1:end-1) - ...
                     (c(1:end-1,2:end)+d(1:end-1,2:end))./C2(1:end-1,2:end);
g = .5*h*g;