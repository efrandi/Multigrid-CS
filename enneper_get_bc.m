function [bottom,top,left,right] = enneper_get_bc(n)
% Computes boundary conditions for the Enneper minimum surface problem.
% INPUT:  n --> the number of grid intevals
% OUTPUT: bottom, top, left, right --> (n+1)-dimensional arrays containing 
%         the boundary values
%
% For each grid point (x,y), the boundary condition is bv(x,y) = u^2 - v^2,
% where u, v are the unique solutions of the equations
%
%  u + u*v^2 - u^3/3 = 0
% -v - u^2*v + v^3/3 = 0

maxit = 5;   % Max iterations in Newton's method
tol = 1e-10; % Tolerance for Newton's method

% Grid parameters
h = 1/n;
x0 = -0.5; xend = 0.5;
y0 = -0.5; yend = 0.5;

% Preallocate the output vectors
bottom = zeros(n+1,1); top = zeros(n+1,1);
left = zeros(n+1,1); right = zeros(n+1,1);
u = [0;0]; nf = [0;0];

for j = 1:4
    % Set up the initial point
    if j == 1 || j == 3, 
        x = x0; y = y0;
    elseif j == 2 
        x = x0; y = yend;
    elseif j == 4 
        x = xend; y = y0;
    end

    % Use Newton's method in R^2 to solve the equations on each grid point
    for i = 1:n+1
        u(1) = x; u(2) = -y; % Initial guess
        for k = 1:maxit
            % Current function value
            nf(1) = u(1) + u(1)*u(2)^2 - u(1)^3/3 - x;
            nf(2) = -u(2) - u(1)^2*u(2) +u(2)^3/3 - y;
            % Stopping criterion on the norm of the nonlinear function
            if norm(nf) <= tol, break; end
            % Define the Jacobian matrix and perform the Newton iteration
            njac(1,1) = 1 + u(2)^2 - u(1)^2;
            njac(1,2) = 2*u(1)*u(2);
            njac(2,1) = -2*u(1)*u(2);
            njac(2,2) = -1 - u(1)^2 + u(2)^2;
            u = u - njac\nf;
        end
        
        if j == 1,
            bottom(i) = u(1)^2 - u(2)^2;
            x = x + h; % Move x to the next grid point, y is fixed
        elseif j == 2
            top(i) = u(1)^2 - u(2)^2;
            x = x + h; % Move x to the next grid point, y is fixed
        elseif j == 3
            left(i) = u(1)^2 - u(2)^2;
            y = y + h; % Move y to the next grid point, x is fixed
        elseif j == 4
            right(i) = u(1)^2 - u(2)^2;
            y = y + h; % Move x to the next grid point, y is fixed
        end
    end
end
    
end