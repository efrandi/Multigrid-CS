function g = morebv_g(x)

n = size(x,1)-1;
h = 1/n;
e = ones(n-1,1);
A = spdiags([-e 2*e -e], -1:1, n-1, n-1);
Ax = x(2:end-1,2:end-1)*A+A*x(2:end-1,2:end-1); %Prodotti di Kronecker
AAx = Ax*A+A*Ax;

tt = linspace(h,1-h,n-1);
[T,S] = meshgrid(tt,tt);
B = x(2:end-1,2:end-1)+T+S+1;
B3 = B.^3;
AB3 = B3*A+A*B3;
B2Ax = B.^2.*Ax;
g = zeros(n+1,n+1);
g(2:end-1,2:end-1) = 2*AAx + h^2*AB3 + 1.5*h^4*B.^5 + 3*h^2*B2Ax;