function g = dpjb_g(x)
%Gradiente obiettivo - problema DPJB

n = size(x,1)-1;
hx = 2*pi/n; hy = 20/n;
xx = linspace(0,2*pi,n+1)';
epsil = 0.1;
g = zeros(n+1,n+1);
wq = repmat((1+epsil*cos(xx)).^3,1,n+1);
wl = repmat(epsil*sin(xx),1,n+1);
al = hx*hy/6*(wq(1:end-1,1:end-1)+wq(2:end,1:end-1)+wq(1:end-1,2:end));
au = hx*hy/6*(wq(2:end,2:end)+wq(1:end-1,2:end)+wq(2:end,1:end-1));
a = x(2:end-1,2:end-1)-x(2:end-1,3:end); b = x(2:end-1,2:end-1)-x(3:end,2:end-1);
c = x(2:end-1,2:end-1)-x(2:end-1,1:end-2); d = x(2:end-1,2:end-1)-x(1:end-2,2:end-1);
g(2:end-1,2:end-1) = al(2:end,2:end).*(b/hx^2+a/hy^2) + au(1:end-1,1:end-1).*(d/hx^2+c/hy^2) + ...
    al(1:end-1,2:end).*d/hx^2 + au(1:end-1,2:end).*a/hy^2 + ...
    al(2:end,1:end-1).*c/hy^2 + au(2:end,1:end-1).*b/hx^2 - hx*hy*wl(2:end-1,2:end-1);