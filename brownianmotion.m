function W=brownianmotion(N,tspan)

T=numel(tspan);

% Simulation of brownianmotion noise that is gaussian distributed
% increments in the concentration of glucose.

randn('state',100) % set the state of randn

dt = T/N; 
dW = zeros(1,N);        % preallocate arrays ...
W = zeros(1,N);         % for efficiency
dW(1) = sqrt(dt)*randn; % first approximation outside the loop ...
W(1) = dW(1);           % since W(0) = 0 is not allowed

    for j = 2:N
    dW(j) = sqrt(dt)*randn; % general increment
    W(j) = W(j-1) + dW(j); % Glucose niveauet stiger gang hvis man plusser med noget her.
    end
       
end
