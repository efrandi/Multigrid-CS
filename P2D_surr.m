function incr = P2D_surr(z,xt,i,j,n,deltaP,deltaM)

incr = [0;0];
h = 1/n;
t = (i-1)*h; s = (j-1)*h;
bij = 2*((1-6*s^2)*t^2*(1-t^2)+(1-6*t^2)*s^2*(1-s^2));

xt(1:3,1:3) = xt(1:3,1:3)-z(1:3,1:3);
%positive search direction
in1 = n^2*deltaP*(-z(1,2)-z(3,2)+4*z(2,2)-z(2,1)-z(2,3)+2*deltaP)-deltaP*bij;
in2 = -n^2*deltaP*(-xt(1,2)-xt(3,2)+4*xt(2,2)-xt(2,1)-xt(2,3)-2*deltaP)+deltaP*bij;
incr(1) = 0.5*(in1+in2);
%negative search direction
in1 = -n^2*deltaM*(-z(1,2)-z(3,2)+4*z(2,2)-z(2,1)-z(2,3)-2*deltaM)+deltaM*bij;
in2 = n^2*deltaM*(-xt(1,2)-xt(3,2)+4*xt(2,2)-xt(2,1)-xt(2,3)+2*deltaM)-deltaM*bij;
incr(2) = 0.5*(in1+in2);

% %positive search direction
% in1 = deltaP*(-z(i-1,j)-z(i+1,j)+4*z(i,j)-z(i,j-1)-z(i,j+1)+2*deltaP)-h^2*deltaP*bij;
% in2 = -deltaP*(-xt(i-1,j)+z(i-1,j)-xt(i+1,j)+z(i+1,j)+4*xt(i,j)-4*z(i,j)...
%       -xt(i,j-1)+z(i,j-1)-xt(i,j+1)+z(i,j+1)-2*deltaP)+h^2*deltaP*bij;
% incr(1) = 0.5*(in1+in2);
% %negative search direction
% in1 = -deltaM*(-z(i-1,j)-z(i+1,j)+4*z(i,j)-z(i,j-1)-z(i,j+1)-2*deltaM)+h^2*deltaM*bij;
% in2 = deltaM*(-xt(i-1,j)+z(i-1,j)-xt(i+1,j)+z(i+1,j)+4*xt(i,j)-4*z(i,j)...
%       -xt(i,j-1)+z(i,j-1)-xt(i,j+1)+z(i,j+1)+2*deltaM)-h^2*deltaM*bij;
% incr(2) = 0.5*(in1+in2);