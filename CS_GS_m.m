function [z,delta,np,fz,ind,kout,grad,unsuccCount] = ...
         CS_GS_m(z,f,update,delta,toll,kmax,rid,aug,xpand,xxtilde)

%CS_GS Gauss-Seidel Coordinate Search with gradient approximation
%
%USAGE:
%   [z,delta,np,fz,ind,kout] = CS_GS(z,f,update,delta,toll,kmax,rid,aug,xpand)
%INPUT:
%   z: initial guess as a (n+1)x(n+1) matrix
%   f: function handle for f
%   update: coordinate-wise fast update for f
%   delta: initial stepsize
%   toll: stopping criterion tolerance
%   kmax: max number of outer iterations
%   rid, aug: reduction/expansion factors for the stepsize
%   xpand: flag for step expansion strategy
%OUTPUT:
%   z: last computed iterate
%   delta: final stepsize
%   np: number of visited trial points
%   fz: final objective function value
%   ind = 0 --> Success (stopping criterion satisfied)
%   ind = -1 --> Failure (kmax reached)
%   kout: number of performed iterations
%   grad: approximation of the gradient

global problem

ind = 0; np = 0; gamma = 0e-4; fz = f(z.sol); 
n = size(z.sol,1)-1; N = (n-1)^2; 
theta = 2; %Expansion parameter
grad = zeros(n+1,n+1);
unsuccCount = 0;

for kout = 1 : kmax,
    %z_old = z.sol;
    succ = 0;
    while succ == 0
        
        if nargin < 10 %No xtilde
            for ik = 1:N %cycle through coordinates 1,...,N
                ix = mod(ik-1,n-1)+2;
                iy = floor((ik-1)/(n-1))+2;
                distP = z.ubound(ix,iy)-z.sol(ix,iy);
                distM = z.sol(ix,iy)-z.lbound(ix,iy);
                stepP = min(delta,distP);
                stepM = min(delta,distM);
%                 ztrial = z.sol;
%                 ztrial(ix,iy) = z.sol(ix,iy) + stepP;
%                 inctrueP = f(ztrial)-f(z.sol);
%                 ztrial(ix,iy) = z.sol(ix,iy) - stepM;
%                 inctrueM = f(ztrial)-f(z.sol);
                
                if strcmp(problem,'morebv') || strcmp(problem,'tutorial');
                    z_window = zeros(5,5);
                    z_window(2:4,2:4) = z.sol(ix-1:ix+1,iy-1:iy+1);
                    if ix < n, z_window(5,2:4) = z.sol(ix+2,iy-1:iy+1); end
                    if ix > 2, z_window(1,2:4) = z.sol(ix-2,iy-1:iy+1); end
                    if iy < n, z_window(2:4,5) = z.sol(ix-1:ix+1,iy+2); end    
                    if iy > 2, z_window(2:4,1) = z.sol(ix-1:ix+1,iy-2); end
                else
                    z_window = z.sol(ix-1:ix+1,iy-1:iy+1);
                end
                incr = update(z_window,ix,iy,n,stepP,stepM);
                
%                 errp = abs(inctrueP-incr(1)); errm = abs(inctrueM-incr(2));
%                 if errp > 1e-13 
%                     fprintf('ErrP = %.4e, (i,j) = (%d,%d) \n',errp,ix,iy);  
%                 end
%                 if errm > 1e-13
%                     fprintf('ErrM = %.4e, (i,j) = (%d,%d) \n',errm,ix,iy);  
%                 end
                
                grad(ix,iy) = (incr(1) - incr(2)) / (stepP+stepM);
                np = np+2;

                [incrmin,imin] = min([incr(1) incr(2)]);
                if incrmin < -gamma*delta^2
                    succ = 1;
                    if imin == 1, step = stepP; dist = distP; else step = -stepM; dist = distM; end
                    if xpand == 1
                        incrnew = update(z_window,ix,iy,n,theta*step,-theta*step);
                        np = np+2;
                        while incrmin < -gamma*step^2 && incrnew(1) < min(-gamma*(theta*step)^2,incrmin) && theta*abs(step) < dist;
                            incrmin = incrnew(1); step = theta*step;
                            incrnew = update(z_window,ix,iy,n,theta*step,-theta*step);
                            np = np+2;
                        end
                    end
                    z.sol(ix,iy) = z.sol(ix,iy)+step; fz = fz+incrmin;
                end
            end
            
        else
            
            for ik = 1:N %cycle through coordinates 1,...,N
                ix = mod(ik-1,n-1)+2;
                iy = floor((ik-1)/(n-1))+2;
                distP = z.ubound(ix,iy)-z.sol(ix,iy);
                distM = z.sol(ix,iy)-z.lbound(ix,iy);
                stepP = min(delta,distP);
                stepM = min(delta,distM);
                if strcmp(problem,'morebv') || strcmp(problem,'tutorial');
                    z_window = zeros(5,5); xt_window = zeros(5,5);
                    z_window(2:4,2:4) = z.sol(ix-1:ix+1,iy-1:iy+1);
                    xt_window(2:4,2:4) = xxtilde(ix-1:ix+1,iy-1:iy+1);
                    if ix < n, z_window(5,2:4) = z.sol(ix+2,iy-1:iy+1);
                               xt_window(5,2:4) = xxtilde(ix+2,iy-1:iy+1); end
                    if ix > 2, z_window(1,2:4) = z.sol(ix-2,iy-1:iy+1); 
                               xt_window(1,2:4) = xxtilde(ix-2,iy-1:iy+1); end
                    if iy < n, z_window(2:4,5) = z.sol(ix-1:ix+1,iy+2); 
                               xt_window(2:4,5) = xxtilde(ix-1:ix+1,iy+2); end    
                    if iy > 2, z_window(2:4,1) = z.sol(ix-1:ix+1,iy-2);
                               xt_window(2:4,1) = xxtilde(ix-1:ix+1,iy-2); end
                else
                    z_window = z.sol(ix-1:ix+1,iy-1:iy+1);
                    xt_window = xxtilde(ix-1:ix+1,iy-1:iy+1);
                end
                incr = update(z_window,xt_window,ix,iy,n,stepP,stepM);
                grad(ix,iy) = (incr(1) - incr(2)) / (stepP+stepM);
                np = np+2;

                [incrmin,imin] = min([incr(1) incr(2)]);
                if incrmin < -gamma*delta^2
                    succ = 1;
                    if imin == 1, step = stepP; dist = distP; else step = -stepM; dist = distM; end
                    if xpand == 1
                        incrnew = update(z_window,xt_window,ix,iy,n,theta*step,-theta*step);
                        np = np+2;
                        while incrmin < -gamma*step^2 && incrnew(1) < min(-gamma*(theta*step)^2,incrmin) && theta*abs(step) < dist;
                            incrmin = incrnew(1); step = theta*step;
                            incrnew = update(z_window,xt_window,ix,iy,n,theta*step,-theta*step);
                            np = np+2;
                        end
                    end
                    z.sol(ix,iy) = z.sol(ix,iy)+step; fz = fz+incrmin;
                end
            end
            
        end
            
        if succ == 0,
            fz = f(z.sol); %Refresh
            delta = rid*delta; unsuccCount = unsuccCount+1;
            if delta < toll, return; end
        else
            delta = aug*delta;
        end    
        
    end
end
%if nargin < 10, fprintf('||x_nu-x_nu-1||_inf = %.4e \n',norm(z.sol(:)-z_old(:),inf)); end
ind = -1;