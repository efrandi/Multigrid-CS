function [f,g] = tutorial_obj(x)

lambda = 10;
n = sqrt(length(x))+1; h = 1/n;
x = reshape(x,n-1,n-1);
e = ones(n-1,1);
A = spdiags([-e 2*e -e], -1:1, n-1, n-1);
Ax = x*A+A*x; %Prodotti di Kronecker
AAx = Ax*A+A*Ax;

tt = linspace(h,1-h,n-1);
[T,S] = meshgrid(tt,tt);
B = lambda*x.*exp(x);
RHS = ((9*pi^2+lambda*exp((T.^2-T.^3).*sin(3*pi*S))).*(T.^2-T.^3)+6*T-2).*sin(3*pi*S);
f = norm(Ax + h^2*(B-RHS),'fro')^2;

C = h^2*(B-RHS);
D = h^2*(B+lambda*exp(x));
AC = C*A+A*C;

g = 2*(AAx + AC + D.*Ax + D.*C); g = g(:);