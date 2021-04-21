function incr = P2D_update(z,i,j,n,deltaP,deltaM)

incr = [0;0];
h = 1/n;
t = (i-1)*h; s = (j-1)*h;
bij = 2*((1-6*s^2)*t^2*(1-t^2)+(1-6*t^2)*s^2*(1-s^2));

%positive search direction
incr(1) = n^2*deltaP*(-z(1,2)-z(3,2)+4*z(2,2)-z(2,1)-z(2,3)+2*deltaP)-deltaP*bij;
%negative search direction
incr(2) = -n^2*deltaM*(-z(1,2)-z(3,2)+4*z(2,2)-z(2,1)-z(2,3)-2*deltaM)+deltaM*bij;

% %positive search direction
% incr(1) = deltaP*(-z(i-1,j)-z(i+1,j)+4*z(i,j)-z(i,j-1)-z(i,j+1)+2*deltaP)-h^2*deltaP*bij;
% %negative search direction
% incr(2) = -deltaM*(-z(i-1,j)-z(i+1,j)+4*z(i,j)-z(i,j-1)-z(i,j+1)-2*deltaM)+h^2*deltaM*bij;