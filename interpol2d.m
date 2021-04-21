function v = interpol2d(u)

%INTERPOL2D Operatore di interpolazione lineare 2-D per il multigrid.
%
%UTILIZZO:
%   v = interpol2d(u)
%   v --> vettore (n+1)^2
%   u --> vettore (nH+1)^2
%   dove nH = n/2.

nH = size(u,1)-1;
n = 2*nH;
v = zeros(n+1,n+1);
v(3:2:n-1,3:2:n-1) = u(2:nH,2:nH);
v(2:2:n,3:2:n-1) = .5*(u(1:nH,2:nH)+u(2:nH+1,2:nH));
v(3:2:n-1,2:2:n) = .5*(u(2:nH,1:nH)+u(2:nH,2:nH+1));
v(2:2:n,2:2:n) = .25*(u(1:nH,1:nH)+u(2:nH+1,1:nH)+u(1:nH,2:nH+1)+u(2:nH+1,2:nH+1));