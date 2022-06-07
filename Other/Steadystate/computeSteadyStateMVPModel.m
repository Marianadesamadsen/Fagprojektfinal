function  [xs, us, flag] = computeSteadyStateMVPModel(t, p, Gs)
% 
% computeSteadyStateMVPModel()
% 
% DESCRIPTION:
% This function solves the equations from the mvpModelSteadyStateWrapper to
% find the steady state vector and the manipulated inputs at steady state. 
% It uses the function fsolve in matlab to solve the euqations when sat to 0.
% To find an appropriate initial guess it uses the naturalMVPmodel.
% 
% INPUT:
%   t  - time
%   p  - a vector parameters               
%   Gs - the steady state blood glucose concentration 
% 
% OUTPUT:
%   xs   - the steady state vector              
%   us   - the steady state manipulated inputs 
%   flag - the flag (integer) returned by fsolve
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

% fsolve options
opts = optimoptions('fsolve',       ...
    'OptimalityTolerance',  1e-10,  ...
    'Display',              'off');

% The disturbance variables
d = 0; % No meal assumed

% Initial guess of the states
xs0 = naturalMVPmodel(p); 

% Initial guess of the manipulated inputs
us0 = [0; 0];

% Initial guess
ws0 = [xs0; us0];

% Solve for the steady state vector and manipulated inputs at steady state
[ws, ~, flag] = fsolve(@mvpModelSteadyStateWrapper, ...
    ws0, opts, t, d, p, Gs);


% OUTPUT

% Steady state vector
xs = ws(1:end-2);

% Steady state manipulated inputs
us = ws(end-1:end);

end


