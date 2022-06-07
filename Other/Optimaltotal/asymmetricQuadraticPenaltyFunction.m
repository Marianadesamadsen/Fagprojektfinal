function phi = asymmetricQuadraticPenaltyFunction(T, Z)
% asymmetricQuadraticPenaltyFunction() 
%
% DESCRIPTION:
% The function evaluates the objective function, which is the integral of
% the asymmetric quadratic glucose penalty function. This integral is 
% computed by using the right-rectangle rule on the given time interval.
%
% INPUT
%   T   - vector of times       (dimension: N+1)
%   Z   - matrix of outputs     (dimension: nz x N+1)
% 
% OUTPUT
%   phi - value of the objective function

% Optimal blood glucose concentration
Gbar = 108.0; % [mg/dL]

% Soft lower bound for the glucose concentration
Gmin = 70; % [mg/dL]

% Hypoglycemia constant
k = 1.0e6;

% Blood glucose concentration to vector form from a matrix
G = Z(:);

% Calculating the penalty function values
rho = 0.5*(G - Gbar).^2 + 0.5*k*max(0.0, (Gmin - G)).^2;

% Approximate the integral using a right-rectangle rule given diff
phi = sum(rho(2:end).*diff(T(:)));

end


