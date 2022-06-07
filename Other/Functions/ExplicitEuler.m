function [T,X] = ExplicitEuler(f, tspank, x0, u, d, p)
%
% ExplicitEuler
%
% DESCRIPTION:
% This function uses the explicit Euler method to approximate a solution
% to a system of differential equations to the initial value problem.
%
% INPUT:
%   f      - function to handle that computes the derivative (MVP model).
%   tspank - points in time where the solution is approximated.
%   x0     - Initial state vector.
%   u      - manipulated inputs (used for MVP model).
%   d      - disturbance inputs (used for MVP model).
%   p      - parameter vector   (used for MVP model).
%
% OUTPUT:
%   T - Time steps of control intervals
%   X - The solution to the differential equations at each control step
% 
% PROJECT:
% Fagprojekt 2022
% A diabetes case study - Meal detection
%
% GENEREL:
% BSc                       : Mathematics and technology
% University                : The Technical University of Denmark (DTU)
% Department                : Applied Mathematics and Computer Science
%
% AUTHORS:
% Emma Victoria Lind
% Mariana de Sá Madsen
% Mona Saleem
%
% CONTACT INFORMATION
% s201205@student.dtu.dk
% s191159@student.dtu.dk
% s204226@student.dtu.dk
%

% Number of control steps
Nk = numel(tspank) - 1;

% Length of the vector
nx = numel(x0);

% Allocating memory for Outputs
T = zeros(Nk+1, 1);
X = zeros(Nk+1, nx);

% Storing initial condition
T(1) = tspank(1);
X(1,:) = x0;

% Step size in the approximation (Explicit Euler step)
h = (tspank(end)-tspank(1))/Nk;

% Overwriting such that we start with tk
tk = tspank(1);

% Overwriting such that we start with xk
xk = x0;

for k=1:Nk

    % Calculating fk (finding derivative with MVP model)
    fk = feval(f, tk, xk, u, d, p);

    % Calculating xk+1 (Explicit Euler step)
    xkp1 = xk + h*fk;

    % Storing it in the matrix
    X(k+1,:) = xkp1;

    % Updating the xk
    xk=xkp1;

    % Updating to tkp1
    tkp1 = tspank(k+1);

    % Storing in vector
    T(k+1) = tkp1;

end

end
