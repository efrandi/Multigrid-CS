function g = finDiff(z,delta)
%FD for P2D

n = size(z.sol,1)-1; N = (n-1)^2;
g = zeros(n+1,n+1);

for ik = 1:N %cycle through coordinates 1,...,N
    ix = mod(ik-1,n-1)+2;
    iy = floor((ik-1)/(n-1))+2;
    stepP = delta; stepM = delta;
    incr = update(z.sol,ix,iy,n,stepP,stepM);
    g(ix,iy) = (incr(1) - incr(2)) / (stepP+stepM);
end

function incr = update(z,i,j,n,deltaP,deltaM)

incr = [0;0];
h = 1/n;
t = (i-1)*h; s = (j-1)*h;
bij = 2*((1-6*s^2)*t^2*(1-t^2)+(1-6*t^2)*s^2*(1-s^2));

%positive search direction
incr(1) = n^2*deltaP*(-z(i-1,j)-z(i+1,j)+4*z(i,j)-z(i,j-1)-z(i,j+1)+2*deltaP)-deltaP*bij;
%negative search direction
incr(2) = -n^2*deltaM*(-z(i-1,j)-z(i+1,j)+4*z(i,j)-z(i,j-1)-z(i,j+1)-2*deltaM)+deltaM*bij;
end

end
