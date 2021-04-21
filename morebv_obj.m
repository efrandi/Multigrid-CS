function [f,g] = morebv_obj(x)

n = sqrt(length(x))+1; h = 1/n;
x = reshape(x,n-1,n-1);

e = ones(n-1,1);
L = spdiags([-e 2*e -e], -1:1, n-1, n-1);
Ax = x*L+L*x; %Prodotti di Kronecker

tt = linspace(h,1-h,n-1);
[T,S] = meshgrid(tt,tt);
B = x+T+S+1;
f = norm(Ax + 0.5*h^2*B.^3,'fro')^2;

AAx = Ax*L+L*Ax;
B3 = B.^3;
AB3 = B3*L+L*B3;
B2Ax = B.^2.*Ax;
g = 2*AAx + h^2*AB3 + 1.5*h^4*B.^5 + 3*h^2*B2Ax; g = g(:);