function [xh,cost,delta] = FMGOpt_m(xh,f,g,delta,Params,update,surr)

%FMGOPT Schema FML/CS
%
%UTILIZZO:
%       [xh,cost,delta] = FMGOpt(xh,f,g,delta,Params)
%INPUT:  
%       xh = soluzione iniziale livello fine
%       f = funzione obiettivo 
%       g = gradiente obiettivo
%       delta = delta iniziale per il CS
%       Params.nu0 = numero v-cicli da fare ad ogni chiamata
%       Params.nu1 = numero cicli di pre-smoothing
%       Params.nu2 = numero cicli di post-smoothing
%       Params.tol = tolleranza per il CS
%       Params.maxit = max cicli di iterazioni per il CS
%       Params.red, Params.aug = parametri di aumento/riduzione del delta nel CS
%       Params.fmgstep = flag per l'utilizzo dell'euristica sui delta
%       update = funzione update obiettivo lungo le coordinate
%OUTPUT:
%       xh = soluzione finale livello fine
%       cost = totale valutazioni di funzione 
%       delta = delta finale in uscita dai v-cicli

n = size(xh.sol,1)-1; %Numero di intervalli interni lungo una sola direzione
l = round(log2(n)); %Livello di discretizzazione, n = 2^l
firstR = xh.sol(1,:); lastR = xh.sol(end,:); 
firstC = xh.sol(:,1); lastC = xh.sol(:,end);


if l <= Params.lmin + 1e-12 %Se siamo sul livello meno fine, soluzione "esatta"
    t0_FM = tic;
    %Proiezione sull'insieme ammissibile
    uViol = find(xh.sol > xh.ubound); lViol = find(xh.sol < xh.lbound);
    xh.sol(uViol) = xh.ubound(uViol); xh.sol(lViol) = xh.lbound(lViol);
    
    % Definizione smoother
    if strcmp(Params.smoother,'J') == 1
        [xh,delta,cost] = CS_J_m(xh,f,update,delta,Params.tol,Params.maxit,Params.red,Params.aug,Params.expand);
    elseif strcmp(Params.smoother,'GS') == 1;
        [xh,delta,cost] = CS_GS_m(xh,f,update,delta,Params.tol,Params.maxit,Params.red,Params.aug,Params.expand);
    else
        error('Please input a valid CS smoother, ''J'' = Jacobi or ''GS'' = Gauss-Seidel');
    end
    t_FM = toc(t0_FM);
    %fprintf('Level %d total cycle time: %.2e \n',l,t_FM);
    return
end

%Altrimenti scendiamo di griglia in maniera ricorsiva

xH.sol = restrict2d(xh.sol);
xH.ubound = xh.ubound(1:2:end,1:2:end); xH.lbound = xh.lbound(1:2:end,1:2:end);

[xH,costH,deltaH] = FMGOpt_m(xH,f,g,delta,Params,update,surr);

cost = costH;

xh.sol = interpol2d(xH.sol); %Le componenti al bordo vanno assegnate
xh.sol(1,:) = firstR; xh.sol(end,:) = lastR; 
xh.sol(:,1) = firstC; xh.sol(:,end) = lastC;

%Proiezione sull'insieme ammissibile
uViol = find(xh.sol > xh.ubound); lViol = find(xh.sol < xh.lbound);
xh.sol(uViol) = xh.ubound(uViol); xh.sol(lViol) = xh.lbound(lViol);

%Definizione euristica FMG
switch Params.fmgstep
    case 0 %Nessuna
        deltanew = 1;
    case 1 %Strategia vecchia (Creta 2012)
        if l == Params.lmin + 1, deltanew = 1; else 
        deltanew = max(Params.tol,min(1,Params.C*deltaH)); 
        end
    case 2 %Strategia nuova (Lisbona 2013)
        deltanew = Params.tol/(Params.C^(l-Params.lmin));
    otherwise 
        error('Please choose a valid strategy');
end

t0_FM = tic;
for i = 1:Params.nu0
    [xh,costV,delta] = VOpt_m(xh,f,g,deltanew,Params,update,surr);
    cost = cost + costV;
end
t_FM = toc(t0_FM);
%fprintf('Level %d total cycle time: %.2e \n',l,t_FM);