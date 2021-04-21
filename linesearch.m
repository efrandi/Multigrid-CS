function [xNew,alpha,fevals,ind] = linesearch(x0,f,d)

%LINESEARCH Line search - decrescita sufficiente con backtracking e espansione.
%
%UTILIZZO:
%   [x,alpha,fevals,ind] = linesearch(x,f,d)
%INPUT: 
%   x0 = punto di partenza, f = obiettivo, d = direzione di ricerca.
%OUTPUT: 
%   xNew = punto finale, alpha = passo finale, fevals = # valutazioni di
%   funzione, ind = 0 successo, -1 fallimento.

ind = 0;
fevals = 0;
alpha = 1; 
alphamax = 64;
rho = 0.5; eta = 2; %Parametri di riduzione e espansione
fx = f(x0);
c = sign(fx)*1e-12; 

xNew = x0 + alpha*d;

%Backtracking (se non ho decrescita provo alpha = rho*alpha, rho^2*alpha...)
while f(xNew) > fx*(1-c) && alpha > eps
    alpha = rho*alpha;
    xNew = x0 + alpha*d;
    fevals = fevals + 1;
end

if alpha <= eps
   ind = -1; xNew = x0;
end

%Espansione (se ho decrescita provo alpha = eta*alpha, eta^2*alpha,...)
if alpha == 1
    fx = f(xNew);
    while f(x0+eta*alpha*d) <= fx*(1-c) && alpha <= alphamax
        alpha = eta*alpha;
        xNew = x0 + alpha*d; fx = f(xNew);
        fevals = fevals + 1;
    end 
end


