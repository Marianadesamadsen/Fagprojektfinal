function [xs, us] = naturalMVPmodel(p)
%
% naturalMVPmodel()
% 
% DESCRIPTION:
% Computes the natural steadystate of the MVP model
% corresponding to no meals or insulin or glucagon. Used as an initial
% guess in computeSteadyStateMVPModel.
%
% INPUT:
% p - vector with parameters                  
%
% OUTPUT:
% xs - the steady states                      
% us - the steady state manipulated inputs    
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

% Parameters
GEZI = p(6);        % [1/min]          Glucose effectiveness
EGP0 = p(7);        % [(mg/dL)/min]    Endogenous glucose production

% Meal subsystem
D1  = 0;            % [g CHO]
D2  = 0;            % [g CHO]

% Insulin subsystem
Isc = 0;            % [mU/L]           Subcutaneous insulin concentration
Ip  = 0;            % [mU/L]           Plasma insulin concentration

% Glucose subsystem
Ieff = 0;           % [1/min]          Insulin effect
G    = EGP0/GEZI;   % [mg/dL]          Blood glucose concentration

% Glucose concentration measurement
Gsc = G;            % [mg/dL]          Subcutaneous glucose concentration

% Basal and bolus insulin flow rates
uba = 0;            % [ mU/min]         Insulin basal flow rate
ubo = 0;            % [ mU/min]         Insulin bolus flow rate

%OUTPUT

% Steady state vector
xs = [D1; D2; Isc; Ip; Ieff; G; Gsc];

% Steady state manipulated inputs
us = [uba; ubo];

end

