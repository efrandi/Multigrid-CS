function f = morebv_f(x)

n = size(x,1)-1;
h = 1/n;
e = ones(n-1,1);
L = spdiags([-e 2*e -e], -1:1, n-1, n-1);
Ax = x(2:end-1,2:end-1)*L+L*x(2:end-1,2:end-1); %Prodotti di Kronecker

tt = linspace(h,1-h,n-1);
[T,S] = meshgrid(tt,tt);
B = x(2:end-1,2:end-1)+T+S+1;
f = norm(Ax + 0.5*h^2*B.^3,'fro')^2;