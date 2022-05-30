function W=brownianmotion(N,T)

randn('state',100) % set the state of randn


dt = T/N; 
dW = zeros(1,N);        % preallocate arrays ...
W = zeros(1,N);         % for efficiency
dW(1) = sqrt(dt)*randn; % first approximation outside the loop ...
W(1) = dW(1);           % since W(0) = 0 is not allowed

    for j = 2:N
    dW(j) = sqrt(dt)*randn; % general increment
    W(j) = W(j-1) + dW(j);
    end

end
