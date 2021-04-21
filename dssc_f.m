function f = dssc_f(x,lambda)
%Funzione obiettivo - problema DSSC

n = size(x,1)-1; h2 = 1/n^2;
fquad = sum(sum((x(2:end-1,2:end)-x(2:end-1,1:end-1)).^2)) + ...
        sum(sum((x(1:end-1,2:end-1)-x(2:end,2:end-1)).^2));
fexp = sum(sum(exp(x(2:end-1,2:end-1)))) + 4*(n-1) + 6/3;
f = .5*fquad-h2*lambda*fexp;
