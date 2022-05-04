function value=singleShootingObjective(ubo,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK)

% Meal and meal bolus 
U(idxbo, 1) = ubo; 

% Simulate 
[T, X] = OpenLoopSimulation(x0, tspan, U, D, p, simModel, simMethod, NK);

% Evaluate the outputs 
Z = CGMsensor(X, p); 

% Evaluate the objective function 
value = abs(scalingFactor*objectiveFunction(T, Z)); 

end
