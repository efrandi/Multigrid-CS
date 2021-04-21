n = 2^4; %Numero intervalli di griglia in ciascuna dimensione
Params.smoother = 'J'; %J = CS-Jacobi, GS = CS-Gauss-Seidel
Params.expand = 1; %Flag espansione nel CS: 1 --> espansione, altrimenti --> no espansione
Params.fmgstep = 2;

%Parametri CS
delta = 1;
Params.red = .25; Params.aug = 1; Params.tol = 1e-5; Params.maxit = 1000;
%Parametri multigrid
Params.nu1 = 2; Params.nu2 = 2; 
Params.lmin = 3;

ncycles = 10;

fprintf('\n');
fprintf('Problem: P2D, Size: %d, Solver: %s, Tol: %.4e \n',...
         (n-1)^2,Params.smoother,Params.tol);

xx = linspace(0,1,n+1);
[X,Y] = meshgrid(xx,xx);

x0.sol = zeros(n+1,n+1); x0.ubound = zeros(n+1,n+1); x0.lbound = zeros(n+1,n+1);
%x0.sol(2:end-1,2:end-1) = ones(n-1,n-1); %Soluzione iniziale
xsin = 0.5*(sin(2*pi*X)+sin(16*pi*X));
ysin = 0.5*(sin(2*pi*Y)+sin(16*pi*Y));
x0.sol = xsin.*ysin;

%Definisco obiettivo, gradiente, update, condizioni al bordo, C euristica

e = ones(n-1,1);
L = spdiags([-e 2*e -e], -1:1, n-1, n-1);
A = n^2*(kron(L,eye(n-1))+kron(eye(n-1),L));
%Soluzione esatta
b = 2*((1-6*X.^2).*Y.^2.*(1-Y.^2)+(1-6*Y.^2).*X.^2.*(1-X.^2)); b = b(2:end-1,2:end-1);
h = 1/n; rhs = 1*b(:);
xes = zeros(n+1,n+1); xes(2:end-1,2:end-1) = reshape(A\rhs,n-1,n-1);
xAnalit = (X.^2-X.^4).*(Y.^4-Y.^2); 
errDiscInf = norm(xes(:)-xAnalit(:),'inf');
errDiscL2D = h*norm(xes(:)-xAnalit(:));
fprintf('Discretization error: inf-norm %.4e, L2D-norm %.4e \n',errDiscInf,errDiscL2D);

x_int = [0 1]; y_int = [0 1];
x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
f = @(x) P2D_f(x); g = @(x) P2D_g(x);
update = @(z,i,j,dim,stepP,stepM) P2D_update(z,i,j,dim,stepP,stepM);
surr = @(z,xt,i,j,dim,stepP,stepM) P2D_surr(z,xt,i,j,dim,stepP,stepM);

Params.approx = 1;
xF = x0; err = h*norm(xes(:)-xF.sol(:));
for j = 1:ncycles
    [xF,costTot] = VOpt_m(xF,f,g,delta,Params,update,surr);
    err_old = err;
    err = h*norm(xes(:)-xF.sol(:));
    fprintf('\n'); fprintf('Cycle %d \n',j);
    fprintf('Fval: %.4e, L2D-error: %.4e, Ratio: %.2e \n',f(xF.sol),err,err/err_old);
end
fprintf('\n');