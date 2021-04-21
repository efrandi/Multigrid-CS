function incr = morebv_surr(z,xt,i,j,n,deltaP,deltaM)

incr = [0;0];
h = 1/n; h2 = h^2;

Az1 = 4*z(3,3)-z(2,3)-z(4,3)-z(3,2)-z(3,4);
c1 = 1+(i-1)*h+(j-1)*h+z(3,3); a1 = Az1+0.5*h2*c1^3;
b = 4*deltaP+0.5*h2*(deltaP^3+3*c1^2*deltaP+3*c1*deltaP^2);
incr11 = b*(b+2*a1);
b = -4*deltaM+0.5*h2*(-deltaM^3-3*c1^2*deltaM+3*c1*deltaM^2);
incr21 = b*(b+2*a1);

if i < n, Az2 = 4*z(4,3)-z(3,3)-z(5,3)-z(4,2)-z(4,4);
    c2 = 1+i*h+(j-1)*h+z(4,3); a2 = Az2+0.5*h2*c2^3;
    incr11 = incr11+deltaP*(deltaP-2*a2); 
    incr21 = incr21+deltaM*(deltaM+2*a2);
end
if i > 2, Az3 = 4*z(2,3)-z(1,3)-z(3,3)-z(2,2)-z(2,4);
    c3 = 1+(i-2)*h+(j-1)*h+z(2,3); a3 = Az3+0.5*h2*c3^3;
    incr11 = incr11+deltaP*(deltaP-2*a3); 
    incr21 = incr21+deltaM*(deltaM+2*a3);
end
if j < n, Az4 = 4*z(3,4)-z(2,4)-z(4,4)-z(3,3)-z(3,5);
    c4 = 1+(i-1)*h+j*h+z(3,4); a4 = Az4+0.5*h2*c4^3;
    incr11 = incr11+deltaP*(deltaP-2*a4); 
    incr21 = incr21+deltaM*(deltaM+2*a4);
end
if j > 2, Az5 = 4*z(3,2)-z(2,2)-z(4,2)-z(3,1)-z(3,3);
    c5 = 1+(i-1)*h+(j-2)*h+z(3,2); a5 = Az5+0.5*h2*c5^3;
    incr11 = incr11+deltaP*(deltaP-2*a5); 
    incr21 = incr21+deltaM*(deltaM+2*a5);
end

xt(2:4,2:4) = xt(2:4,2:4)-z(2:4,2:4);

Az1 = 4*xt(3,3)-xt(2,3)-xt(4,3)-xt(3,2)-xt(3,4);
c1 = 1+(i-1)*h+(j-1)*h+xt(3,3); a1 = Az1+0.5*h2*c1^3;
b = -4*deltaP+0.5*h2*(-deltaP^3-3*c1^2*deltaP+3*c1*deltaP^2);
incr12 = b*(b+2*a1);
b = 4*deltaM+0.5*h2*(deltaM^3+3*c1^2*deltaM+3*c1*deltaM^2);
incr22 = b*(b+2*a1);

if i < n, xt(5,3) = xt(5,3)-z(5,3); 
    Az2 = 4*xt(4,3)-xt(3,3)-xt(5,3)-xt(4,2)-xt(4,4);
    c2 = 1+i*h+(j-1)*h+xt(4,3); a2 = Az2+0.5*h2*c2^3;
    incr12 = incr12+deltaP*(deltaP+2*a2); 
    incr22 = incr22+deltaM*(deltaM-2*a2);
end
if i > 2, xt(1,3) = xt(1,3)-z(1,3);
    Az3 = 4*xt(2,3)-xt(1,3)-xt(3,3)-xt(2,2)-xt(2,4);
    c3 = 1+(i-2)*h+(j-1)*h+xt(2,3); a3 = Az3+0.5*h2*c3^3;
    incr12 = incr12+deltaP*(deltaP+2*a3); 
    incr22 = incr22+deltaM*(deltaM-2*a3);
end
if j < n, xt(3,5) = xt(3,5)-z(3,5);
    Az4 = 4*xt(3,4)-xt(2,4)-xt(4,4)-xt(3,3)-xt(3,5);
    c4 = 1+(i-1)*h+j*h+xt(3,4); a4 = Az4+0.5*h2*c4^3;
    incr12 = incr12+deltaP*(deltaP+2*a4); 
    incr22 = incr22+deltaM*(deltaM-2*a4);
end
if j > 2, xt(3,1) = xt(3,1)-z(3,1);
    Az5 = 4*xt(3,2)-xt(2,2)-xt(4,2)-xt(3,1)-xt(3,3);
    c5 = 1+(i-1)*h+(j-2)*h+xt(3,2); a5 = Az5+0.5*h2*c5^3;
    incr12 = incr12+deltaP*(deltaP+2*a5); 
    incr22 = incr22+deltaM*(deltaM-2*a5);
end

incr(1) = 0.5*(incr11+incr12);
incr(2) = 0.5*(incr21+incr22);