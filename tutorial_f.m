function f = tutorial_f(x)

lambda = 10;
n = size(x,1)-1;
h = 1/n;
e = ones(n-1,1);
L = spdiags([-e 2*e -e], -1:1, n-1, n-1);
Ax = x(2:end-1,2:end-1)*L+L*x(2:end-1,2:end-1); %Prodotti di Kronecker

tt = linspace(h,1-h,n-1);
[T,S] = meshgrid(tt,tt);
B = lambda*x(2:end-1,2:end-1).*exp(x(2:end-1,2:end-1));
RHS = ((9*pi^2+lambda*exp((T.^2-T.^3).*sin(3*pi*S))).*(T.^2-T.^3)+6*T-2).*sin(3*pi*S);
f = norm(Ax + h^2*(B-RHS),'fro')^2;