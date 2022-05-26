function [T,X] = OpenLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, NK)
%
% OPENLOOPSIMULATION Open-loop simulation
%
% DESCRIPTION:
% Perform an open-loop simulation for given initial condition, control
% intervals, disturbance variables, parameters, and simulation model.
%
% INPUT:
%   x0                  - initial state                                     (dimension: nx    )
%   tspan               - boundaries of the control intervals               (dimension: N+1   )
%   U                   - disturbance variables for each control interval   (dimension: nu x N)
%   D                   - disturbance variables for each control interval   (dimension: nd x N)
%   p                   - parameters                                        (dimension: np    )
%   simModel            - simulation model          (function handle)
%   simMethod           - simulation method         (function handle)
%   Nk                  - Number of time steps in each control interval 
% 
% OUTPUT:
%   T - The control state of time for each stept                            (dimension:      N+1)
%   X - The states of the solved differential equation stores in matrix     (dimension: nx x N+1)
%

% Number of control steps 
N = numel(tspan) - 1; 

% Number of states 
nx = numel(x0); 

% Number of time steps in each control interval 
Nk=NK;

% Allocate memory 
T = zeros(1, N+1); 
X = zeros(nx, N+1); 

% Initial condition in each control interval
xk = x0; 

% Store solution
T(1) = tspan(1);
X(:,1) = x0;

for k = 1:N 
    % Times 
    tk    = tspan(k);
    tkp1  = tspan(k+1);
    
    % Manipulated inputs and disturbance variables
    uk = U(:,k); 
    dk = D(:,k); 
    
    % Time interval
    tspank = linspace(tk, tkp1, Nk+1);
    
    % Solve initial value problem 
    [Tk, Xk] = simMethod(simModel, tspank, xk, uk, dk, p);
    
    % Update initial condition 
    xk = Xk(end, :)'; 
    
    % Store solution 
    T(k+1) = Tk(end)';
    X(:, k+1) = Xk(end, :)';
    
end

     
end

