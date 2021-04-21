function f = mins_f(x)

n = size(x,1)-1;
h = 1/n;
a = x(1:end-1,1:end-1)-x(1:end-1,2:end);
b = x(2:end,2:end)-x(1:end-1,2:end);
c = x(2:end,2:end)-x(2:end,1:end-1);
d = x(1:end-1,1:end-1)-x(2:end,1:end-1);
C1 = sqrt(h^2 + a.^2 + b.^2);
C2 = sqrt(h^2 + c.^2 + d.^2);
f = .5*h*sum(sum(C1+C2));