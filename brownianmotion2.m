function W=brownianmotion2(Nk,tspan)

tk=tspan(end)-tspan(1);

%randn('state',100) % set the state of randn

dt = tk/Nk;
dW = sqrt(dt)*randn(1,Nk); % increments
W = cumsum(dW); % cumulative sum

end