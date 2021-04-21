function g = tutorial_g(x)

lambda = 10;
n = size(x,1)-1;
h = 1/n;
e = ones(n-1,1);
A = spdiags([-e 2*e -e], -1:1, n-1, n-1);
Ax = x(2:end-1,2:end-1)*A+A*x(2:end-1,2:end-1); %Prodotti di Kronecker
AAx = Ax*A+A*Ax;

tt = linspace(h,1-h,n-1);
[T,S] = meshgrid(tt,tt);
B = lambda*x(2:end-1,2:end-1).*exp(x(2:end-1,2:end-1));
RHS = ((9*pi^2+lambda*exp((T.^2-T.^3).*sin(3*pi*S))).*(T.^2-T.^3)+6*T-2).*sin(3*pi*S);
C = h^2*(B-RHS);
D = h^2*(B+lambda*exp(x(2:end-1,2:end-1)));
AC = C*A+A*C;

g = zeros(n+1,n+1);
g(2:end-1,2:end-1) = 2*(AAx + AC + D.*Ax + D.*C);