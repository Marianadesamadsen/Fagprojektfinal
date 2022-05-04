function f = MVPmodel(t, x, u, d, p)
% MVPMODEL Evaluate the right-hand side function of the Medtronic Virtual
% Patient (MPV) model.
% 
% Description:
% This model describes the dynamics of the blood glucose concentration in
% response to meals and subcutaneous insulin.
%
% Auther: 
% Mariana de Sá Madsen 
% Mona Saleem 
% Emma Lind
%
% INPUT:
%   t - time
%   x - a vector of state variables         (dimension:  7)
%   u - a vector of manipulated inputs      (dimension:  2)
%   d - a vector of disturbance variables   (dimension:  1)
%   p - a vector parameters                 (dimension: 10)
% 
% OUTPUT:
% Returns all the devivatives of X.


% Inizializing the values of the vector X.

% Meal subsystem
D1  = x(1); % [g CHO]
D2  = x(2); % [g CHO]

% Insulin subsystem
Isc = x(3); % [mU/L] Subcutaneous insulin concentration (Under the skin)
Ip  = x(4); % [mU/L] Plasma insulin concentration (In the blood)
Ieff = x(5); % [1/min]  Effekt of the insulin

% Glucose subsystem
G = x(6); % [mg/dL] Blood glucose concentration

% Glucose concentration measurement
Gsc = x(7); % [mg/dL] Subcutaneous glucose concentration (Under the skin)

% Manipulated inputs
ubasal = u(1); % [ mU/min] Insulin basal flow rate
ubolus = u(2); % [ mU/min] Insulin bolus flow rate

% Parameters
tau1    = p(1); % [min]            Insulin absorption time
tau2    = p(2); % [min]            Insulin absorption time
CI      = p(3); % [dL/min]         Insulin clearance rate
p2      = p(4); % [1/min]          Inverse of insulin action time constant
SI      = p(5); % [(dL/mU)/min]    Insulin sensitivity
GEZI    = p(6); % [1/min]          Glucose effectiveness
EGP0    = p(7); % [(mg/dL)/min]    Endogenous glucose production
VG      = p(8); % [dL]             Glucose distribution volume
taum    = p(9); % [min]            Meal time constant
tausc   = p(10); % [min]           Subcutaneous insulin time constant

% Meal rate of appearance
RA = 1000*D2/(VG*taum); % [(mg/dL)/min]

% Inizialising the output vector
f = zeros(7, 1);

% Calculating the derivatives:
f(1) =  d  - D1 /taum;
f(2) = (D1 - D2)/taum;
f(3) = ((ubasal + ubolus)/CI - Isc)/tau1;
f(4) = (Isc-Ip)/tau2;
f(5) = p2*(SI*Ip - Ieff);
f(6) = -(GEZI + Ieff)*G + EGP0 + RA;
f(7) = (G - Gsc)/tausc;

end



