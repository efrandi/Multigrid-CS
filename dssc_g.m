function g = dssc_g(x,lambda)
%Gradiente obiettivo - problema DSSC

n = size(x,1)-1; h2 = 1/n^2;
g = zeros(n+1,n+1);
gquad = 4*x(2:end-1,2:end-1)-x(1:end-2,2:end-1)-x(3:end,2:end-1)-x(2:end-1,1:end-2)-x(2:end-1,3:end);
gexp = exp(x(2:end-1,2:end-1));
g(2:end-1,2:end-1) = gquad-h2*lambda*gexp;