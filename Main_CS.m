%MAIN - ALGORITMO CS

n = 2^5; %Numero intervalli di griglia in ciascuna dimensione
Params.smoother = 'J'; %J = CS-Jacobi, GS = CS-Gauss-Seidel
Params.expand = 0; %Flag espansione nel CS: 1 --> espansione, altrimenti --> no espansione
global problem
problem = 'p2d';

%Parametri CS
delta = 1;
Params.red = .25; Params.aug = 1; Params.tol = 1e-4; 
Params.maxit = 50000;

fprintf('\n');
fprintf('Problem: %s, Size: %d, Solver: %s, Tol: %.4e \n',...
         problem,(n-1)^2,Params.smoother,Params.tol);

x0.sol = zeros(n+1,n+1); x0.ubound = zeros(n+1,n+1); x0.lbound = zeros(n+1,n+1);
x0.sol(2:end-1,2:end-1) = ones(n-1,n-1); %Soluzione iniziale

xx = linspace(0,1,n+1);
[X,Y] = meshgrid(xx,xx);

%Definisco obiettivo, gradiente, update, condizioni al bordo, C euristica
switch problem
    case 'p2d'
        
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
        f = @(x) P2D_f(x); g = @(x) P2D_g(x); Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) P2D_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) P2D_surr(z,xt,i,j,dim,stepP,stepM);
                
    case 'mins-sb'
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) mins_f(x); g = @(x) mins_g(x); Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) mins_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) mins_surr(z,xt,i,j,dim,stepP,stepM);
        bc = xx.*(1-xx); x0.sol(1,:) = bc; x0.sol(end,:) = bc;
       
    case 'dssc'
        lambda = 5;
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) dssc_f(x,lambda); g = @(x) dssc_g(x,lambda); Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) dssc_update(z,i,j,dim,stepP,stepM,lambda);
        surr = @(z,xt,i,j,dim,stepP,stepM) dssc_surr(z,xt,i,j,dim,stepP,stepM,lambda);
    
    case 'morebv'
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) morebv_f(x); g = @(x) morebv_g(x); Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) morebv_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) morebv_surr(z,xt,i,j,dim,stepP,stepM);
        
    case 'dpjb'
        x_int = [0 2*pi]; y_int = [0 20];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = 0;
        f = @(x) dpjb_f(x); g = @(x) dpjb_g(x); Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) dpjb_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) dpjb_surr(z,xt,i,j,dim,stepP,stepM);
    
    otherwise
        error('Please input a valid problem name');
end


tmain = tic;
if strcmp(Params.smoother,'J') == 1
    [xF,delta,evals,fval,~,iters] = CS_J_m(x0,f,update,delta,Params.tol,Params.maxit,Params.red,Params.aug,Params.expand);
elseif strcmp(Params.smoother,'GS') == 1
    [xF,delta,evals,fval,~,iters] = CS_GS_m(x0,f,update,delta,Params.tol,Params.maxit,Params.red,Params.aug,Params.expand);
else 
    error('Please define a valid CS strategy, ''J'' = Jacobi or ''GS'' = Gauss-Seidel');
end
time1 = toc(tmain);
fprintf('\n'); fprintf('Fevals: %.4e Iters: %d Fval: %.15e Time: %.2e \n',evals,iters,fval,time1);
if strcmp(problem,'p2d');
    fprintf('Inf-error: %.4e, L2D-error: %.4e \n',norm(xes(:)-xF.sol(:),'inf'),h*norm(xes(:)-xF.sol(:)));
end