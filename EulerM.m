function [T,X] = EulerM(f,tspan,x0,u,d,p)

% Number of control steps
N = numel(tspan) - 1;

% Length of the vector
nx = numel(x0);

% Allocating memory for Outputs
T = zeros(N+1, 1);
X = zeros(N+1, nx);

% Storing initial condition
T(1) = tspan(1);
X(1,:) = x0;

% Overwriting such that we start with tk
tk = tspan(1);

% Overwriting such that we start with xk
xk = x0;

% The stochastic elements 
lambda = 1.5;
mu = 0; 

% 
R = 4; 

for k=1:N

    dt = T(k)/(2^8);
    Dt = R*dt; 
    
    dW = sqrt(dt)*randn(1,2^8);
    
    % Calculating fk (finding derivative with MVP model)
    fk = feval(f, tk, xk, u, d, p);

    Winc = sum(dW(R*(k-1)+1:R*k));
    Xtemp = fk + Dt*lambda*fk + mu*fk*Winc;
    X(k) = Xtemp;

    % Storing it in the matrix
    X(k+1,:) = xkp1;

    % Updating the xk
    xk=xkp1;

    % Updating to tkp1
    tkp1 = tspan(k+1);

    % Storing in vector
    T(k+1) = tkp1;
   
end


end


