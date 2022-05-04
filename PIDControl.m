function [uk,ctrlstate] = PIDControl(yk, ctrlPar, ctrlState)
%
% PIDCONTROLLER PID controller for controlling the insulin flow rate.
%
% DESCRIPTION: 
% This function implements a discretized proportional-integral-derivative
% (PID) controller for controlling the blood glucose concentration.
%
% INPUT:
%   Ik              - the integral term (Ik)
%   ybar            - the traget glucose concentration, y=108
%   yk              - the current glucose concentration, yk
%   ykm1            - previous glucose concentration, yk-1
%   us              - insulin state
%   Kp,Ti,Td        - tuned parameters
%   Ts              - samlping time, 5 min
% 
% Output:
%   uk          - a vector of manipulated inputs
%   crtlstate   - the updated controller state


% Unpack control parameters
Ts      = ctrlPar(1); % [min]    Sampling time
Kp      = ctrlPar(2); %          Proportional gain
% Ki      = ctrlPar(3); %          Integrator gain
% Kd      = ctrlPar(4); %          Derivative gain
ybar    = ctrlPar(5); % [mg/dL]  Target blood glucose concentration
ubar    = ctrlPar(6); % [mU/min] Nominal insulin flow rate
Ti      = ctrlPar(7);
Td      = ctrlPar(8);

% Unpack control state
Ik = ctrlState(1); %           Value of integral at previous time step
ykm1 = ctrlState(2); % [mg/dL]   Previous observed glucose concentration


% Computing

ek = yk-ybar; % Setpoint error

Ki = Kp * Ts/Ti; % helps controlling the steady state
Kd = Kp * Td/Ts; % the top 

Pk = Kp * ek; % Proportional term: controls how fast we change the error 

Ikp1 = Ik + Ki * ek; % Integral term: Ik the area of the error

Dk = Kd * (yk-ykm1); % Derivative term: The top of the curve

uba = ubar + Pk + Ik + Dk; % Basal insulin flow rate

ubo = 0; % Bolus insulin flow rate is zero

% OUTPUT 

% Manipulated inputs (must be non-negative) OUTPUT
uk = [uba,ubo];

% Controller state OUTPUT
ctrlstate = [Ikp1; yk];

end








