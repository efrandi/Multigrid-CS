function u = enneper_get_sol(n)
% Computes the analytic solution of the Enneper minimum surface problem.
% INPUT:  n --> the number of grid intevals
% OUTPUT: u --> (n+1)^2-dimensional array containing the solution
%
% For each grid point (x,y), the solution is u(x,y) = v(1)^2 - v(2)^2,
% where v the unique solutions of the equations
%
%  v(1) + v(1)*v(2)^2 - v(1)^3/3 = 0
% -v(2) - v(1)^2*v(2) + v(2)^3/3 = 0

maxit = 10;   % Max iterations in Newton's method
tol = 1e-12; % Tolerance for Newton's method

% Grid parameters
h = 1/n;
x0 = -0.5; x = x0;
y0 = -0.5; y = y0;

% Preallocate the solution
u = zeros(n+1,n+1);
v = [0;0]; nf = [0;0];

for i = 1:n+1
    for j = 1:n+1
        % Use Newton's method to solve the equations on each grid point
        v(1) = x; v(2) = -y; % Initial guess
        for k = 1:maxit
            % Current function value
            nf(1) = v(1) + v(1)*v(2)^2 - v(1)^3/3 - x;
            nf(2) = -v(2) - v(1)^2*v(2) +v(2)^3/3 - y;
            % Stopping criterion on the norm of the nonlinear function
            if norm(nf) <= tol, break; end
            % Define the Jacobian matrix and perform the Newton iteration
            njac(1,1) = 1 + v(2)^2 - v(1)^2;
            njac(1,2) = 2*v(1)*v(2);
            njac(2,1) = -2*v(1)*v(2);
            njac(2,2) = -1 - v(1)^2 + v(2)^2;
            v = v - njac\nf;
        end
        
        u(i,j) = v(1)^2 - v(2)^2;
        x = x + h; % Move x to the next grid point, y is fixed
    end
    x = x0;    %Restart x
    y = y + h; % Move y to the next grid point, y is fixed
end
    
end