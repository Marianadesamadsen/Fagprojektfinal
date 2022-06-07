
clear

%% Loading all folders
% 
fprintf('Loading diabetes library .. ');

% Add real thermodynamics functions
addpath(genpath(fullfile(pwd, './Other')));

% Let the user know that the library is being loaded
fprintf('Done\n');

%% Script start

tspank = linspace(0, 5, 10+1);
Nk = numel(tspank) - 1;

W=brownianmotion(Nk,tspank);

figure 
plot(W,'*')

% Den viser 100 målinger på et sekund. Den viser random gå tur rundt
% omkring 0. 

% W er den random process i målingerne. 


