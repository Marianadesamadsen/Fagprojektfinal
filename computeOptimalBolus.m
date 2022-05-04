function ubo=computeOptimalBolus(N,M,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK)

% INPUT
% N: Ubo  
% M: Interval
% The others are for singleshootingobjective

% Calculating the max phi value 
min_phi=singleShootingObjective(N,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK);

% Setting optimal bolus for that value
ubo=N; 

for i=1:M:N
    
    % Calculating new objective function values given the interval of ubo
    % in the loop
current_phi=singleShootingObjective(i,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK);

% Checking if the current value is less than the min_phi
if abs(current_phi)<abs(min_phi)   
    ubo=i; % Setting ubo to that value
    min_phi=current_phi; % Overwriting for testing new current_phi
end

end



