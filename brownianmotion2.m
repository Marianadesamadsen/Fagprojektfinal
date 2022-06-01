function W=brownianmotion2(N,tspan)

T=tspan(end);

randn('state',100) % set the state of randn

dt = T/N;
dW = sqrt(dt)*randn(1,N); % increments
W = cumsum(dW); % cumulative sum

end