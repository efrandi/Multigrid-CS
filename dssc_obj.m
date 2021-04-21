function [f,g] = dssc_obj(x_vec)
%Funzione obiettivo - problema DSSC
lambda = 5;

n = sqrt(length(x_vec))+1; h2 = 1/n^2;
x = zeros(n+1,n+1);
x(2:end-1,2:end-1) = reshape(x_vec,n-1,n-1);

fquad = sum(sum((x(2:end-1,2:end)-x(2:end-1,1:end-1)).^2)) + ...
        sum(sum((x(1:end-1,2:end-1)-x(2:end,2:end-1)).^2));
fexp = sum(sum(exp(x(2:end-1,2:end-1)))) + 4*(n-1) + 6/3;
f = .5*fquad-h2*lambda*fexp;

gquad = 4*x(2:end-1,2:end-1)-x(1:end-2,2:end-1)-x(3:end,2:end-1)-x(2:end-1,1:end-2)-x(2:end-1,3:end);
gexp = exp(x(2:end-1,2:end-1));
g = gquad-h2*lambda*gexp; g = g(:);