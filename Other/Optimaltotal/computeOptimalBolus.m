function ubo=computeOptimalBolus(N,M,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK)
% computeOptimalBolus()
% 
% DESCRIPTION:
% COMPUTEOPTIMALBOLUS compute an optimal insulin bolus for a given
% objective function and simulation model.
% This function uses fmincon to solve a dynamic optimization problem where
% the decision variable is the insulin bolus administered during the first
% control interval. The decision variable is chosen such as to minimize a
% user-provided objective function over a user-defined prediction horizon.
% The objective function, the simulation model, the output model, the
% initial condition, the disturbances, and the model parameters are
% provided by the user.
%
% INPUT:
% N                 - ubo, bolus insulin                                (scalar)
% M                 - Interval
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
% singleShootingObjective()
% 
% OUTPUT:
% ubo               - the local optimal bolus insulin 
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

% Calculating the min phi value 
min_phi=singleShootingObjective(N,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK);

% Setting optimal bolus for the highest posible value
ubo=N; 

for i=1:M:N
    
    % Calculating new objective function values given the interval of ubo
    % in the loop
current_phi=singleShootingObjective(i,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK);

% Checking if the current value is less than the min_phi
if abs(current_phi)<abs(min_phi)   
    ubo=i; % Setting ubo to that value that is lower than ubo before
    min_phi=current_phi; % Overwriting for testing new current_phi
end

end



