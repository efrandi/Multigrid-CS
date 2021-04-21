function v = restrict2d(u)

%RESTRICT2D Operatore di restrizione "full weighting" 2-D per il multigrid.
%
%UTILIZZO:
%   v = restrict2d(u)
%   u --> vettore (n+1)^2
%   v --> vettore (nH+1)^2
%   dove nH = n/2.

n = size(u,1)-1;
nH = n/2;
v = zeros(nH+1,nH+1);
v(2:end-1,2:end-1) = (1/16)*(4*u(3:2:n-1,3:2:n-1) + ...
2*(u(3:2:n-1,2:2:n-2)+u(3:2:n-1,4:2:n)+u(2:2:n-2,3:2:n-1)+u(4:2:n,3:2:n-1)) ...
+u(2:2:n-2,2:2:n-2)+u(2:2:n-2,4:2:n)+u(4:2:n,2:2:n-2)+u(4:2:n,4:2:n));

%Assegno le condizioni al bordo
v(1,:) = u(1,1:2:end);
v(end,:) = u(end,1:2:end);
v(:,1) = u(1:2:end,1); 
v(:,end) = u(1:2:end,end);