
clear
tspank = linspace(0, 5, 10+1);
Nk = numel(tspank) - 1;

W=brownianmotion(Nk,tspank);

figure 
plot(W)

% Den viser 100 målinger på et sekund. Den viser random gå tur rundt
% omkring 0. 

% W er den random process i målingerne. 


