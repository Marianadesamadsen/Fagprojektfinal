function [xs, us] = naturalMVPmodel(p)
%
% DESCRIPTION: Computes the natrual steadystate of the MVP model
% corresponding to no meals or insulin or glucagon derived analytically 
%
% INPUT: 
% p - vector with parameters                  (dim(p)=31)
% 
% OUTPUT:
% xs - the steady states                      (dim(xs)=16)
% us - the steady state manipulated inputs    (dim(us)=3)
%

% Parameters
GEZI = p(6); % [1/min]          Glucose effectiveness
EGP0 = p(7); % [(mg/dL)/min]    Endogenous glucose production

% Meal subsystem
D1  = 0; % [g CHO]
D2  = 0; % [g CHO]

% Insulin subsystem
Isc = 0; % [mU/L] Subcutaneous insulin concentration
Ip  = 0; % [mU/L] Plasma insulin concentration

% Glucose subsystem
Ieff = 0;         % [1/min] Insulin effect
G    = EGP0/GEZI; % [mg/dL] Blood glucose concentration

% Glucose concentration measurement
Gsc = G; % [mg/dL] Subcutaneous glucose concentration

% Basal and bolus insulin flow rates
uba = 0; % [ mU/min] Insulin basal flow rate
ubo = 0; % [ mU/min] Insulin bolus flow rate

% Steady state
xs = [D1; D2; Isc; Ip; Ieff; G; Gsc];

% Steady state manipulated input
us = [uba; ubo];

end

