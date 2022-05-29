function R = mvpModelSteadyStateWrapper(w, t, d, p, Gs)
%
% mvpModelSteadyStateWrapper
% 
% DESCRIPTION:
% The function is used to identidy the steady state. It computes three
% equations. These will be used with a foot-finding algorithm, meaning that
% the equations will be sat equal to 0, and then computing the manipulated
% values and state vector when steady state. 
%
% INPUT:
%   w  - a vector of states and manipulated inputs. [x,u]
%   t  - time
%   d  - a vector of disturbance variables         
%   p  - a vector of parameters                        
%   Gs - the steady state blood glucose concentration [mg/dL]
%
% OUTPUT:
% R - the residual equations 
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

% Extract state vector
x = w(1:end-2);

% Extract manipulated inputs
u = w(end-1:end);

% Glucose value
G = x(6); % [mg/dL]

% Bolus insulin 
ubo = u(2); % [mU/min]

% Computing the right hand side of the MVP model (the derivatives).
R1 = MVPmodel(t, x, u, d, p);

% The change in the blood glucose concentration from the steady state
R2 = G - Gs; % [mg/dL]

% Specification equation for bolus insulin flow rate
R3 = ubo;

% All these equations will be sat equal to 0, and solved with fsolve to 
% compute the state vector and the manipulated inputs at steady state

% Collected residual equations
R = [R1; R2; R3]; 

end