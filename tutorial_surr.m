function incr = tutorial_surr(z,xt,i,j,n,deltaP,deltaM)

lambda = 10;
incr = [0;0];
h = 1/n; h2 = h^2;

y = (i-1)*h; x = (j-1)*h;
f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);

Az1 = 4*z(3,3)-(z(2,3)+z(4,3)+z(3,2)+z(3,4));
nnl = lambda*z(3,3)*exp(z(3,3));
a1 = Az1+h2*(nnl-f);
b = 4*deltaP+h2*(lambda*(z(3,3)+deltaP)*exp(z(3,3)+deltaP)-nnl);
incr11 = b*(b+2*a1);
b = -4*deltaM+h2*(lambda*(z(3,3)-deltaM)*exp(z(3,3)-deltaM)-nnl);
incr21 = b*(b+2*a1);

if i < n, Az2 = 4*z(4,3)-(z(3,3)+z(5,3)+z(4,2)+z(4,4));
    y = i*h; x = (j-1)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(4,3)*exp(z(4,3));
    a2 = Az2+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP-2*a2); 
    incr21 = incr21+deltaM*(deltaM+2*a2);
end
if i > 2, Az3 = 4*z(2,3)-(z(1,3)+z(3,3)+z(2,2)+z(2,4));
    y = (i-2)*h; x = (j-1)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(2,3)*exp(z(2,3));
    a3 = Az3+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP-2*a3); 
    incr21 = incr21+deltaM*(deltaM+2*a3);
end
if j < n, Az4 = 4*z(3,4)-(z(2,4)+z(4,4)+z(3,3)+z(3,5));
    y = (i-1)*h; x = j*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(3,4)*exp(z(3,4));
    a4 = Az4+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP-2*a4); 
    incr21 = incr21+deltaM*(deltaM+2*a4);
end
if j > 2, Az5 = 4*z(3,2)-(z(2,2)+z(4,2)+z(3,1)+z(3,3));
    y = (i-1)*h; x = (j-2)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(3,2)*exp(z(3,2));
    a5 = Az5+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP-2*a5); 
    incr21 = incr21+deltaM*(deltaM+2*a5);
end

xt(2:4,2:4) = xt(2:4,2:4)-z(2:4,2:4);

y = (i-1)*h; x = (j-1)*h;
f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);

Az1 = 4*xt(3,3)-(xt(2,3)+xt(4,3)+xt(3,2)+xt(3,4));
nnl = lambda*xt(3,3)*exp(xt(3,3));
a1 = Az1+h2*(nnl-f);
b = -4*deltaP+h2*(lambda*(z(3,3)-deltaP)*exp(z(3,3)-deltaP)-nnl);
incr12 = b*(b+2*a1);
b = 4*deltaM+h2*(lambda*(z(3,3)+deltaM)*exp(z(3,3)+deltaM)-nnl);
incr22 = b*(b+2*a1);

if i < n, xt(5,3) = xt(5,3)-z(5,3);  
    Az2 = 4*xt(4,3)-(xt(3,3)+xt(5,3)+xt(4,2)+xt(4,4));
    y = i*h; x = (j-1)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*xt(4,3)*exp(xt(4,3));
    a2 = Az2+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP+2*a2); 
    incr21 = incr21+deltaM*(deltaM-2*a2);
end
if i > 2, xt(1,3) = xt(1,3)-z(1,3);
    Az3 = 4*xt(2,3)-(xt(1,3)+xt(3,3)+xt(2,2)+xt(2,4));
    y = (i-2)*h; x = (j-1)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*xt(2,3)*exp(xt(2,3));
    a3 = Az3+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP+2*a3); 
    incr21 = incr21+deltaM*(deltaM-2*a3);
end
if j < n, xt(3,5) = xt(3,5)-z(3,5);
    Az4 = 4*xt(3,4)-(xt(2,4)+xt(4,4)+xt(3,3)+xt(3,5));
    y = (i-1)*h; x = j*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*xt(3,4)*exp(xt(3,4));
    a4 = Az4+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP+2*a4); 
    incr21 = incr21+deltaM*(deltaM-2*a4);
end
if j > 2, xt(3,1) = xt(3,1)-z(3,1);
    Az5 = 4*xt(3,2)-(xt(2,2)+xt(4,2)+xt(3,1)+xt(3,3));
    y = (i-1)*h; x = (j-2)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*xt(3,2)*exp(xt(3,2));
    a5 = Az5+h2*(nnl-f);
    incr11 = incr11+deltaP*(deltaP+2*a5); 
    incr21 = incr21+deltaM*(deltaM-2*a5);
end

incr(1) = 0.5*(incr11+incr12);
incr(2) = 0.5*(incr21+incr22);
