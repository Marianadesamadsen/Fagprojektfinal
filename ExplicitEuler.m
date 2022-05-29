function [T,X] = ExplicitEuler(f, tspan, x0, u, d, p)
%
% ExplicitEUler fixed step size to approximate the solution to 
% an initial value problem.
%
% [T,X] = ExplicitEuler_vores1(f, tspan, x0, varvargin)
% 
% DESCRIPTION:
% Approximating the solution of the differential equations 
% to the initial value problem
%
% INPUT:
%   f      - function handle for evaluating the right-hand side function
%   tspan  - points in time where the solution is approximated (dimension: N+1)
%   x0     - initial states                                    (dimension: nx)
%   u,d,pf - inputs for the function for evaluating 
% 
% OUTPUT:
%   T - boundaries of control intervals             (dimension: N+1)
%   X - The solution to the differential equations  (dimension: N+1 x nx)


% Number of steps
N=numel(tspan) - 1;

% Length of the vector
nx = numel(x0);

% Allocating memory:
T = zeros(N+1, 1);
X = zeros(N+1, nx); 

% Storing initial condition
T(1) = tspan(1);
X(1,:) = x0;

% Step size
h = (tspan(end)-tspan(1))/N;

% Initial time
tk=tspan(1);

% Overwriting such that we start with xk
xk=x0;

for k=1:N
    
    % Calculating fk
    fk = feval(f, tk, xk, u, d, p);
    
    % Calculating xk+1
    xkp1 = xk + h*fk;
    % Storing it in the matrix
    X(k+1,:) = xkp1;
    % Updating the xk
    xk=xkp1; 
    
    % Calculating tkp1
    tkp1 = tspan(k+1);
    % Storing in the vector
    T(k+1) = tkp1;
    
end

end



