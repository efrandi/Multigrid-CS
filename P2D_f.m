function f = P2D_f(v)

n = size(v,1)-1;
h = 1/n;

e = ones(n-1,1);
L = spdiags([-e 2*e -e], -1:1, n-1, n-1);

xx = linspace(h,1-h,n-1);
[X,Y] = meshgrid(xx,xx);
b = zeros(n+1,n+1);
b(2:end-1,2:end-1) = 2*((1-6*X.^2).*Y.^2.*(1-Y.^2)+(1-6*Y.^2).*X.^2.*(1-X.^2));

Av = zeros(n+1,n+1);
%Formula ricavata dai prodotti di Kronecker
Av(2:end-1,2:end-1) = v(2:end-1,2:end-1)*L+L*v(2:end-1,2:end-1); 
f = .5*n^2*sum(sum(Av.*v)) - sum(sum(b.*v));
%f = .5*sum(sum(Av.*v)) - h^2*sum(sum(b.*v)); 