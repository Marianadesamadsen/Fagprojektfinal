function W = brownianmotion(Nk,tspan)
% brownianmotion()
% 
% DESCRIPTION:
% This function simulates a Wiener process. That is making a random walk
% around 0. This function will simulate the noise of a brownian motion that
% will be used for Euler Maruyama. 
%
% INPUT:
% Nk     - Number of time steps in each control/sampling interval
% tspan  - points in time where the solution is approximated.
%
% OUTPUT:
% W      - The brownian motion path
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

tk=(tspan(end)-tspan(1))/Nk;

dt = tk/Nk;
dW = sqrt(dt)*randn(1,Nk); % increments
W = cumsum(dW); % cumulative sum

end
