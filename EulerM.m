function [T,X] = EulerM(f,tspank,x0,u,d,p)

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

% Overwriting such that we start with tk
tk = tspank(1);

% delta t
dt = tspank(2)-tspank(1);

% Overwriting such that we start with xk
xk = x0;

% The simulation of noise
W=brownianmotion2(Nk+1,tspank); % Snak med hjælpelærerne 

for k=1:Nk
    
    % Calculating fk (finding derivative with MVP model)
    fk = feval(f, tk, xk, u, d, p);
    
    gt=zeros(size(xk));
    gt(6)=xk(6);

    % Euler Maruyama step
    xkp1 = xk + fk*dt + gt*(W(k+1)-W(k));

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


