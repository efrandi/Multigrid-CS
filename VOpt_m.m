function [xh,cost,delta] = VOpt_m(xh,f,g,delta,Params,update,surr,vh)

%VOPT Schema ML/CS
%
%UTILIZZO:
%       [xh,cost,delta] = VOpt(xh,f,g,delta,Params,vh)
%INPUT: 
%       xh = soluzione iniziale livello fine
%       f = funzione obiettivo
%       g = gradiente obiettivo
%       delta = delta iniziale
%       Params.nu1 = numero cicli di pre-smoothing
%       Params.nu2 = numero cicli di post-smoothing
%       Params.tol = tolleranza per il CS
%       Params.maxit = max cicli di iterazioni per il CS
%       Params.red, Params.aug = parametri di aumento/riduzione del delta
%       Params.approx = flag per l'utilizzo del surrogato approssimato
%       vh = termine di correzione
%       update = funzione update obiettivo lungo le coordinate
%OUTPUT:
%       xh = soluzione finale livello fine
%       cost = totale valutazioni di funzione
%       delta = delta finale

n = size(xh.sol,1)-1; %Numero di intervalli interni lungo una sola direzione
l = round(log2(n)); %Livello di discretizzazione, n = 2^l

delta0 = delta; %Delta iniziale, usato per far partire il v-ciclo nella ricorsione
tol_solver = delta-1e-12; %Il solver al livello "coarsest" fa una sola riduzione

%Strategia sul criterio di arresto nello smoother
if Params.fmgstep == 2, %Vedo FML/CS come metodo "diretto"
    tol_cs = 0; %Disattivo il criterio di arresto nello smoothing
else %Vedo FML/CS come metodo iterativo
    tol_cs = Params.tol; 
end

%Definizione smoother
if strcmp(Params.smoother,'J') == 1
    CS = @(y1,y2,y3,y4,y5,y6,yopt) CS_J_m(y1,y2,y3,y4,y5,y6,Params.red,Params.aug,Params.expand);
elseif strcmp(Params.smoother,'GS') == 1;
    CS = @(y1,y2,y3,y4,y5,y6,yopt) CS_GS_m(y1,y2,y3,y4,y5,y6,Params.red,Params.aug,Params.expand);
else
    error('Please input a valid CS smoother, ''J'' = Jacobi or ''GS'' = Gauss-Seidel');
end

if nargin < 8 %Chiamata livello piu' fine: niente surrogato, no vincoli
    
    [xh,delta,cost,~,~,~,gh] = CS(xh,f,update,delta,tol_cs,Params.nu1);
    if delta < tol_cs, return; end
%     g_true = g(xh.sol);
%     gh = finDiff(xh,delta);
%     fprintf('Current stepsize = %.4e, ||g_h(xtilde)-gh||_inf = %.4e \n',...
%             delta,norm(g_true(:)-gh(:),inf));
    
    xH.sol = restrict2d(xh.sol);
    %Definizione vincoli x liv. inferiore (Gratton & co. 2011)
    xH.lbound = xH.sol + max(max(xh.lbound(2:end-1,2:end-1)-xh.sol(2:end-1,2:end-1)));
    xH.ubound = xH.sol + min(min(xh.ubound(2:end-1,2:end-1)-xh.sol(2:end-1,2:end-1)));
    if Params.approx == 0, 
        gH = g(xH.sol); rgh = restrict2d(g(xh.sol));
        vH = gH - rgh;
    else
        vH = restrict2d(gh);
    end
    t0_V = tic;
    [xH2,costH] = VOpt_m(xH,f,g,delta0,Params,update,surr,vH);
    t_V = toc(t0_V); %fprintf('    Level %d recursion time: %.2e \n',l-1,t_V);
    cost = cost + costH;
    e = interpol2d(xH2.sol-xH.sol);
    [xh.sol] = linesearch(xh.sol,f,e); 

    [xh,delta,fevals] = CS(xh,f,update,delta,tol_cs,Params.nu2);
    cost = cost + fevals;
    return
end

if Params.approx == 0 %Uso il modello classico
    fs = @(z) f(z) - sum(sum(vh.*z));
    fincr = @(z,i,j,dim,stepP,stepM) update(z,i,j,dim,stepP,stepM) - vh(i,j)*[stepP;-stepM];
    if l <= Params.lmin + 1e-12 %Se siamo sul livello piu' basso, chiamo il solver diretto con tolleranza tol_solver
        [xh,delta,cost] = CS(xh,fs,fincr,delta,tol_solver,Params.maxit);
        if Params.approx ~= 0, cost = 2*cost; end 
        return
    end
    %Altrimenti avvio lo schema di correzione
    [xh,delta,cost] = CS(xh,fs,fincr,delta,tol_cs,Params.nu1);
    if delta < tol_cs, return; end 
    xH.sol = restrict2d(xh.sol);
    xH.lbound = xH.sol + max(max(xh.lbound(2:end-1,2:end-1)-xh.sol(2:end-1,2:end-1)));
    xH.ubound = xH.sol + min(min(xh.ubound(2:end-1,2:end-1)-xh.sol(2:end-1,2:end-1)));
    gH = g(xH.sol); rgh = restrict2d(g(xh.sol));
    vH = restrict2d(vh) + gH - rgh;
    t0_V = tic;
    [xH2,costH] = VOpt_m(xH,f,g,delta0,Params,update,surr,vH);
    %t_V = toc(t0_V); fprintf('    Level %d recursion time: %.2e \n',l-1,t_V);
    cost = cost + costH;
    e = interpol2d(xH2.sol-xH.sol); %Le componenti al bordo nell'interpolazione sono nulle 
    [xh.sol] = linesearch(xh.sol,fs,e);
    [xh,delta,fevals] = CS(xh,fs,fincr,delta,tol_cs,Params.nu2);
    cost = cost + fevals;

else %Uso il modello surrogato approssimato
    
    xxtilde = 2*xh.sol;
    fs = @(z) .5*(f(z)+f(xxtilde-z)) + sum(sum(vh.*z));
    fincr = @(z,xt,i,j,dim,stepP,stepM) surr(z,xt,i,j,dim,stepP,stepM) + vh(i,j)*[stepP;-stepM];
    
    %Ri-definizione smoother
    if strcmp(Params.smoother,'J') == 1
        CS = @(y1,y2,y3,y4,y5,y6) CS_J_m(y1,y2,y3,y4,y5,y6,Params.red,Params.aug,Params.expand,xxtilde);
    elseif strcmp(Params.smoother,'GS') == 1;
        CS = @(y1,y2,y3,y4,y5,y6) CS_GS_m(y1,y2,y3,y4,y5,y6,Params.red,Params.aug,Params.expand,xxtilde);
    else
        error('Please input a valid CS smoother, ''J'' = Jacobi or ''GS'' = Gauss-Seidel');
    end

    if l <= Params.lmin + 1e-12 %Se siamo sul livello piu' basso, chiamo il solver diretto con tolleranza tol_solver
        [xh,delta,cost] = CS(xh,fs,fincr,delta,tol_solver,Params.maxit);
        if Params.approx ~= 0, cost = 2*cost; end 
        return
    end
    %Altrimenti avvio lo schema di correzione
    [xh,delta,cost,~,~,~,gh] = CS(xh,fs,fincr,delta,tol_cs,Params.nu1);
    cost = 2*cost; %Costo doppio del surrogato approx.
    if delta < tol_cs, return; end 
    xH.sol = restrict2d(xh.sol);
    xH.lbound = xH.sol + max(max(xh.lbound(2:end-1,2:end-1)-xh.sol(2:end-1,2:end-1)));
    xH.ubound = xH.sol + min(min(xh.ubound(2:end-1,2:end-1)-xh.sol(2:end-1,2:end-1)));
    vH = restrict2d(gh);
    t0_V = tic;
    [xH2,costH] = VOpt_m(xH,f,g,delta0,Params,update,surr,vH);
    %t_V = toc(t0_V); fprintf('    Level %d recursion time: %.2e \n',l-1,t_V);
    cost = cost + costH;
    e = interpol2d(xH2.sol-xH.sol); %Le componenti al bordo nell'interpolazione sono nulle 
    [xh.sol] = linesearch(xh.sol,fs,e);
    [xh,delta,fevals] = CS(xh,fs,fincr,delta,tol_cs,Params.nu2);
    fevals = 2*fevals;
    cost = cost + fevals;
end