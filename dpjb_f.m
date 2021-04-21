function f = dpjb_f(x)
%Funzione obiettivo - problema DJPB

n = size(x,1)-1;
hx = 2*pi/n; hy = 20/n;
xx = linspace(0,2*pi,n+1)';
epsil = 0.1;
wq = repmat((1+epsil*cos(xx)).^3,1,n+1);
wl = repmat(epsil*sin(xx),1,n+1);
al = hx*hy/6*(wq(1:end-1,1:end-1)+wq(2:end,1:end-1)+wq(1:end-1,2:end));
au = hx*hy/6*(wq(2:end,2:end)+wq(1:end-1,2:end)+wq(2:end,1:end-1));
ql = al.*(1/hx^2*(x(1:end-1,1:end-1)-x(2:end,1:end-1)).^2 + 1/hy^2*(x(1:end-1,1:end-1)-x(1:end-1,2:end)).^2);
qu = au.*(1/hx^2*(x(2:end,2:end)-x(1:end-1,2:end)).^2 + 1/hy^2*(x(2:end,2:end)-x(2:end,1:end-1)).^2);
f = 0.5*(sum(sum(ql))+sum(sum(qu))) - hx*hy*sum(sum(wl.*x));
