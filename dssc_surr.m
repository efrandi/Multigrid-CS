function incr = dssc_surr(z,xt,i,j,n,deltaP,deltaM,lambda)
%Aggiornamento lungo le coordinate - problema DSSC

incr = [0;0];
h = 1/n; 

a = z(2,2)-z(1,2); c = z(2,2)-z(2,1);
b = z(2,2)-z(3,2); d = z(2,2)-z(2,3);
exx = exp(z(2,2));
i11 = 2*deltaP^2 + deltaP*(a+b+c+d) + lambda*h^2*exx*(1-exp(deltaP));
i21 = 2*deltaM^2 - deltaM*(a+b+c+d) + lambda*h^2*exx*(1-exp(-deltaM));

xt(1:3,1:3) = xt(1:3,1:3)-z(1:3,1:3);

a = xt(2,2)-xt(1,2); c = xt(2,2)-xt(2,1);
b = xt(2,2)-xt(3,2); d = xt(2,2)-xt(2,3);
exx = exp(xt(2,2));
i12 = 2*deltaP^2 - deltaP*(a+b+c+d) + lambda*h^2*exx*(1-exp(-deltaP));
i22 = 2*deltaM^2 + deltaM*(a+b+c+d) + lambda*h^2*exx*(1-exp(+deltaM));

incr(1) = 0.5*(i11+i12);
incr(2) = 0.5*(i21+i22);

