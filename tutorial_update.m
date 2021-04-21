function incr = tutorial_update(z,i,j,n,deltaP,deltaM)

lambda = 10;
incr = [0;0];
h = 1/n; h2 = h^2;

y = (i-1)*h; x = (j-1)*h;
f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);

Az1 = 4*z(3,3)-(z(2,3)+z(4,3)+z(3,2)+z(3,4));
nnl = lambda*z(3,3)*exp(z(3,3));
a1 = Az1+h2*(nnl-f);
b = 4*deltaP+h2*(lambda*(z(3,3)+deltaP)*exp(z(3,3)+deltaP)-nnl);
incr(1) = b*(b+2*a1);
b = -4*deltaM+h2*(lambda*(z(3,3)-deltaM)*exp(z(3,3)-deltaM)-nnl);
incr(2) = b*(b+2*a1);

if i < n, Az2 = 4*z(4,3)-(z(3,3)+z(5,3)+z(4,2)+z(4,4));
    y = i*h; x = (j-1)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(4,3)*exp(z(4,3));
    a2 = Az2+h2*(nnl-f);
    incr(1) = incr(1)+deltaP*(deltaP-2*a2); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a2);
end
if i > 2, Az3 = 4*z(2,3)-(z(1,3)+z(3,3)+z(2,2)+z(2,4));
    y = (i-2)*h; x = (j-1)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(2,3)*exp(z(2,3));
    a3 = Az3+h2*(nnl-f);
    incr(1) = incr(1)+deltaP*(deltaP-2*a3); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a3);
end
if j < n, Az4 = 4*z(3,4)-(z(2,4)+z(4,4)+z(3,3)+z(3,5));
    y = (i-1)*h; x = j*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(3,4)*exp(z(3,4));
    a4 = Az4+h2*(nnl-f);
    incr(1) = incr(1)+deltaP*(deltaP-2*a4); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a4);
end
if j > 2, Az5 = 4*z(3,2)-(z(2,2)+z(4,2)+z(3,1)+z(3,3));
    y = (i-1)*h; x = (j-2)*h;
    f = ((9*pi^2+lambda*exp((x^2-x^3)*sin(3*pi*y)))*(x^2-x^3)+6*x-2)*sin(3*pi*y);
    nnl = lambda*z(3,2)*exp(z(3,2));
    a5 = Az5+h2*(nnl-f);
    incr(1) = incr(1)+deltaP*(deltaP-2*a5); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a5);
end
