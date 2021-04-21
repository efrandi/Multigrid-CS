function incr = morebv_update(z,i,j,n,deltaP,deltaM)

incr = [0;0];
h = 1/n; h2 = h^2;

Az1 = 4*z(3,3)-z(2,3)-z(4,3)-z(3,2)-z(3,4);
c1 = 1+(i-1)*h+(j-1)*h+z(3,3); a1 = Az1+0.5*h2*c1^3;
b = 4*deltaP+0.5*h2*(deltaP^3+3*c1^2*deltaP+3*c1*deltaP^2);
incr(1) = b*(b+2*a1);
b = -4*deltaM+0.5*h2*(-deltaM^3-3*c1^2*deltaM+3*c1*deltaM^2);
incr(2) = b*(b+2*a1);

if i < n, Az2 = 4*z(4,3)-z(3,3)-z(5,3)-z(4,2)-z(4,4);
    c2 = 1+i*h+(j-1)*h+z(4,3); a2 = Az2+0.5*h2*c2^3;
    incr(1) = incr(1)+deltaP*(deltaP-2*a2); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a2);
end
if i > 2, Az3 = 4*z(2,3)-z(1,3)-z(3,3)-z(2,2)-z(2,4);
    c3 = 1+(i-2)*h+(j-1)*h+z(2,3); a3 = Az3+0.5*h2*c3^3;
    incr(1) = incr(1)+deltaP*(deltaP-2*a3); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a3);
end
if j < n, Az4 = 4*z(3,4)-z(2,4)-z(4,4)-z(3,3)-z(3,5);
    c4 = 1+(i-1)*h+j*h+z(3,4); a4 = Az4+0.5*h2*c4^3;
    incr(1) = incr(1)+deltaP*(deltaP-2*a4); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a4);
end
if j > 2, Az5 = 4*z(3,2)-z(2,2)-z(4,2)-z(3,1)-z(3,3);
    c5 = 1+(i-1)*h+(j-2)*h+z(3,2); a5 = Az5+0.5*h2*c5^3;
    incr(1) = incr(1)+deltaP*(deltaP-2*a5); 
    incr(2) = incr(2)+deltaM*(deltaM+2*a5);
end