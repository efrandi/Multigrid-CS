function incr = dssc_update(z,i,j,n,deltaP,deltaM,lambda)
%Aggiornamento lungo le coordinate - problema DSSC

incr = [0;0];
h = 1/n; 

a = z(2,2)-z(1,2); c = z(2,2)-z(2,1);
b = z(2,2)-z(3,2); d = z(2,2)-z(2,3);
exx = exp(z(2,2));
incr(1) = 2*deltaP^2 + deltaP*(a+b+c+d) + lambda*h^2*exx*(1-exp(deltaP));
incr(2) = 2*deltaM^2 - deltaM*(a+b+c+d) + lambda*h^2*exx*(1-exp(-deltaM));