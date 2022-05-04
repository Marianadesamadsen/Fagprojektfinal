function  [xs, us, flag] = computeSteadyStateMVPModel(t, p, Gs)
%
% COMPUTESTEADYSTATEMVPMODEL Compute the steady state of the Medtronic
% Virtual Patient (MVP) model.
%
% DESCRIPTION:
% Solve a combination of the nonlinear right-hand side function and a set
% of linear specification equations for the steady state of the Medtronic
% Virtual Patient (MVP) model at given blood glucose concentration and zero
% insulin and glucagon boli:
%
% Here, w = [x; u]. This set of nonlinear and linear equations are solved
% using fsolve to a hardcoded tolerance of 10^(-10).
% 
% INPUT:
%   t  - time
%   p  - a vector parameters               (dimension: 10)
%   Gs - the steady state blood glucose concentration [mg/dL]
%
%
% OUTPUT:
%   xs   - the steady state                    (dimension: 7)
%   us   - the steady state manipulated inputs (dimension: 2)
%   flag - the flag (integer) returned by fsolve


% fsolve options
opts = optimoptions('fsolve',       ...
    'OptimalityTolerance',  1e-10,  ...
    'Display',              'off');

% The disturbance variables
d = 0;

% Initial guess of the states
xs0 = naturalMVPmodel(p);

% Initial guess of the manipulated inputs
us0 = [0; 0];

% Initial guess
ws0 = [xs0; us0];

% Solve for the steady state
[ws, ~, flag] = fsolve(@mvpModelSteadyStateWrapper, ...
    ws0, opts, t, d, p, Gs);

% Steady state
xs = ws(1:end-2);

% Steady state manipulated inputs
us = ws(end-1:end);

end