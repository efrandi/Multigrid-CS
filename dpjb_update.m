function incr = dpjb_update(z,i,j,n,deltaP,deltaM)
%Aggiornamento lungo le coordinate - problema DPJB

incr = [0;0];
hx = 2*pi/n; hy = 20/n; epsil = 0.1;
d = z(i,j)-z(i-1,j); c = z(i,j)-z(i,j-1);
b = z(i,j)-z(i+1,j); a = z(i,j)-z(i,j+1);
wq = (1+epsil*cos((i-1)*hx))^3;
wqu = (1+epsil*cos((i-2)*hx))^3;
wql = (1+epsil*cos(i*hx))^3;
al = hx*hy/6*(2*wq+wql); au = hx*hy/6*(2*wq+wqu);
ali = hx*hy/6*(2*wqu+wq); alj = hx*hy/6*(2*wq+wql);
aui = hx*hy/6*(2*wql+wq); auj = hx*hy/6*(2*wq+wqu);
wlij = epsil*sin((i-1)*hx);
incr(1) = 0.5*deltaP*(al*((deltaP+2*a)/hy^2+(deltaP+2*b)/hx^2)+au*((deltaP+2*c)/hy^2+(deltaP+2*d)/hx^2)+ ...
          ali*(deltaP+2*d)/hx^2+auj*(deltaP+2*a)/hy^2+alj*(deltaP+2*c)/hy^2+aui*(deltaP+2*b)/hx^2) - ...
          hx*hy*wlij*deltaP;
incr(2) = 0.5*deltaM*(al*((deltaM-2*a)/hy^2+(deltaM-2*b)/hx^2)+au*((deltaM-2*c)/hy^2+(deltaM-2*d)/hx^2)+ ...
          ali*(deltaM-2*d)/hx^2+auj*(deltaM-2*a)/hy^2+alj*(deltaM-2*c)/hy^2+aui*(deltaM-2*b)/hx^2) + ...
          hx*hy*wlij*deltaM;