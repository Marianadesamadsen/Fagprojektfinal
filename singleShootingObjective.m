function value = singleShootingObjective(ubo,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK)
% singleShootingObjective()
% 
% DESCRIPTION:
% The function evaluates the objective function for a dynamical
% optimization problem for computing the optimal bolus insulin by using a
% single shooting approach. The model (objective function) is constrained 
% by the initial state vector, the meal and the insulin rates for which the
% the singleshooting function solves for giving the objective function
% value. 
%
% INPUT:
% ubo               - bolus insulin                                     (scalar)
% idxbo             - index for ubo                                     (scalar)
% scalingFactor     - scales the value for the objective function       (scalar)
% x0                - initial state vector                              (dimension 7)
% tspan             - time interval to integrate over                   (dimension N+1)
% U                 - matrix of bolus and basal insulin (manipulated input)
% D                 - matrix of meal rate (disturbance)                 (dimension nd x N)
% p                 - parameter values                                  (dimension np)
% simModel          - simulation model, MVPmodel                        (function handle)
% simMethod         - simulation method, ExplicitEuler                  (function handle)
% objectiveFunction - simulation method/function                        (function handle)
% NK                - Number of timesteps in each time interval         (scalar)
%
% DEPENDENCIES:
% OpenLoopSimulation()
% CGMsensor()
% 
% OUTPUT:
% value - the value for the objective function
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
% Meal and meal bolus 
U(idxbo, 1) = ubo; 

% Simulate 
[T, X] = OpenLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, NK);

% Evaluate the outputs 
Z = CGMsensor(X, p); 

% Evaluate the objective function 
value = abs(scalingFactor*objectiveFunction(T, Z)); 

end




