%MAIN - ALGORITMO FML/CS
%Risolve problemi discretizzati su griglie bidimensionali
%
%PROBLEMI TEST:
%   'p2d': Problema di Poisson in 2D (lineare)
%   'mins-sb': Superficie minima con bordo smooth (nonlineare convessa)
%   'dssc': Modello di combustione solida (quadratica)
%   'morebv': Problema test di More' (non convessa)
%   'dpjb': Pressione lubrificante in una boccola cilindrica (quadratica, vincoli non-negativita')

n = 2^6; %Numero intervalli di griglia in ciascuna dimensione
%Params.approx = 1; %Flag per il surrogato approssimato: 0 --> surr esatto, altrimenti --> surr approx
Params.fmgstep = 2; %Flag euristica sui delta in FML/CS: 0 --> no euristica, 1 --> veloce, 2 --> accurata
Params.smoother = 'GS'; %J = CS-Jacobi, GS = CS-Gauss-Seidel
Params.expand = 1; %Flag espansione nel CS: 1 --> espansione, altrimenti --> no espansione
global problem
problem = 'mins-null';

%Parametri CS
delta = 1;
Params.red = .25; Params.aug = 1; Params.tol = 1e-4; Params.maxit = 10000;
%Parametri multigrid
Params.nu0 = 1; Params.nu1 = 10; Params.nu2 = 10; 
Params.lmin = 3;
Params.C = 8;

fprintf('\n');
fprintf('Problem: %s, Size: %d, Solver: %s, Tol: %.4e \n',...
         problem,(n-1)^2,Params.smoother,Params.tol);
% if Params.approx == 0, fprintf('Classical surrogate \n'); 
% else fprintf('Approximate surrogate \n'); end
%fprintf('\n');

xx = linspace(0,1,n+1);
[X,Y] = meshgrid(xx,xx);

x0.sol = zeros(n+1,n+1); x0.ubound = zeros(n+1,n+1); x0.lbound = zeros(n+1,n+1);

%Soluzione iniziale (1,...,1)'

%x0.sol(2:end-1,2:end-1) = ones(n-1,n-1);

%Combinazione di seni (x esperimenti smoothing)
xsin = 0.5*(sin(2*pi*X)+sin(16*pi*X));
ysin = 0.5*(sin(2*pi*Y)+sin(16*pi*Y));
x0.sol = xsin.*ysin;

%subplot(1,2,1); surf(X,Y,x0.sol); xlim([0 1]); ylim([0 1]);

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
        f = @(x) P2D_f(x); g = @(x) P2D_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) P2D_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) P2D_surr(z,xt,i,j,dim,stepP,stepM);
                
    case 'mins-sb'
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) mins_f(x); g = @(x) mins_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) mins_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) mins_surr(z,xt,i,j,dim,stepP,stepM);
        bc = xx.*(1-xx); x0.sol(1,:) = bc; x0.sol(end,:) = bc;
        
    case 'mins-lin'
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) mins_f(x); g = @(x) mins_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) mins_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) mins_surr(z,xt,i,j,dim,stepP,stepM);
        Z = 1+2*X+3*Y; xes = Z; h = 1/n;
        x0.sol(1,:) = Z(1,:); x0.sol(end,:) = Z(end,:);
        x0.sol(:,1) = Z(:,1); x0.sol(:,end) = Z(:,end);
        
     case 'mins-enneper'
        x_int = [-0.5 0.5]; y_int = [-0.5 0.5];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) mins_f(x); g = @(x) mins_g(x);
        update = @(z,i,j,dim,stepP,stepM) mins_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) mins_surr(z,xt,i,j,dim,stepP,stepM);
        Z = enneper_get_sol(n); xes = Z; h = 1/n;
        x0.sol(1,:) = Z(1,:); x0.sol(end,:) = Z(end,:);
        x0.sol(:,1) = Z(:,1); x0.sol(:,end) = Z(:,end);
        
    case 'mins-null'
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) mins_f(x); g = @(x) mins_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) mins_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) mins_surr(z,xt,i,j,dim,stepP,stepM);
        Z = zeros(n+1,n+1); xes = Z; h = 1/n;
        x0.sol(1,:) = Z(1,:); x0.sol(end,:) = Z(end,:);
        x0.sol(:,1) = Z(:,1); x0.sol(:,end) = Z(:,end);
       
    case 'dssc'
        lambda = 5;
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) dssc_f(x,lambda); g = @(x) dssc_g(x,lambda); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) dssc_update(z,i,j,dim,stepP,stepM,lambda);
        surr = @(z,xt,i,j,dim,stepP,stepM) dssc_surr(z,xt,i,j,dim,stepP,stepM,lambda);
    
    case 'morebv'
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) morebv_f(x); g = @(x) morebv_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) morebv_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) morebv_surr(z,xt,i,j,dim,stepP,stepM);
        
    case 'tutorial'
        xes = (X.^2-X.^3).*sin(3*pi*Y); h = 1/n;
        x_int = [0 1]; y_int = [0 1];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = -inf;
        f = @(x) tutorial_f(x); g = @(x) tutorial_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) tutorial_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) tutorial_surr(z,xt,i,j,dim,stepP,stepM);
        
    case 'dpjb'
        x_int = [0 2*pi]; y_int = [0 20];
        x0.ubound(1:end,1:end) = +inf; x0.lbound(1:end,1:end) = 0;
        f = @(x) dpjb_f(x); g = @(x) dpjb_g(x); %Params.C = 2^2;
        update = @(z,i,j,dim,stepP,stepM) dpjb_update(z,i,j,dim,stepP,stepM);
        surr = @(z,xt,i,j,dim,stepP,stepM) dpjb_surr(z,xt,i,j,dim,stepP,stepM);
    
    otherwise
        error('Please input a valid problem name');
end

% Params.approx = 0; fprintf('\n'); fprintf('Classical surrogate \n'); 
% fprintf('-------------------\n');
% tmain_1 = tic;
% [xF1,costTot] = FMGOpt_m(x0,f,g,delta,Params,update,surr);
% time1 = toc(tmain_1);
% fprintf('\n'); fprintf('Fevals: %.4e Fval: %.15e Time: %.2e \n',costTot,f(xF1.sol),time1);
% if exist('xes','var') == 1
%     fprintf('Inf-error: %.4e, L2D-error: %.4e \n',norm(xes(:)-xF1.sol(:),'inf'),h*norm(xes(:)-xF1.sol(:)));
% end

Params.approx = 1; fprintf('\n'); fprintf('Approximate surrogate \n');
fprintf('---------------------\n');
tmain_2 = tic;
[xF2,costTot] = FMGOpt_m(x0,f,g,delta,Params,update,surr);
time2 = toc(tmain_2);
fprintf('\n'); fprintf('Fevals: %.4e Fval: %.15e Time: %.2e \n',costTot,f(xF2.sol),time2);
if exist('xes','var') == 1    
    fprintf('Inf-error: %.4e, L2D-error: %.4e \n',norm(xes(:)-xF2.sol(:),'inf'),h*norm(xes(:)-xF2.sol(:)));
end
% 
% fprintf('\n'); fprintf('||x_clas - x_appr||_inf = %.4e, t_appr/t_clas = %.2e \n',norm(xF1.sol(:)-xF2.sol(:),'inf'),time2/time1);
fprintf('\n');

% Plot della soluzione
[X,Y] = meshgrid(linspace(x_int(1),x_int(2),n+1),linspace(y_int(1),y_int(2),n+1)); 
surf(X,Y,100*xF2.sol); xlim(x_int); ylim(y_int);