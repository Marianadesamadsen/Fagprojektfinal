
clear all 
clc 
close all

%% Loading all folders
fprintf('Loading diabetes library .. ');

% Add real thermodynamics functions
addpath(genpath(fullfile(pwd, './Other')));

% Let the user know that the library is being loaded
fprintf('Done\n');

%% Formatting the plots 

fs = 11; % Font size 
lw = 3; % Line width 
set(groot, 'DefaultAxesFontSize',   fs); % Set default font size

% Set default line widths  
set(groot, 'DefaultLineLineWidth',  lw);
set(groot, 'DefaultStairLineWidth', lw);
set(groot, 'DefaultStemLineWidth',  lw);

%reset(groot);

%%  Conversion factors 

h2min = 60;      % Convert from h   to min 
min2h = 1/h2min; % Convert from min to h 
U2mU  = 1e3;     % Convert from U   to mU 
mU2U  = 1/U2mU;  % Convert from mU  to U 

%% Inizializing parameters

p = [49,47,20.1,0.0106,0.0081,0.0022,1.33,253,47,5]; % Parameters at steady state
Gs = 108; % [mg/dL]: Steady state blood glucose concentration
ts = []; 

%% Computing steadty state

[xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);

% If fsolve did not converge, throw an error
if(flag ~= 1), error ('fsolve did not converge!'); end


%% Intital and final time for the simulation over 18 hours

t0 =  0;       % min - start time
tf = 18*h2min; % min - end time
Ts = 5;        % min - Sampling time

%% Number of contral intervals

N = (tf - t0)/Ts; % [#]

%% Number of time steps in each control/sampling interval
Nk = 10;

%% Time span
tspan = 5*(0:N);

%% Initial condition

x0 = xs;

% Manipulated inputs
U = repmat(us, 1, N);

% Disturbance variables
D = zeros(1, N);

%% Inizialising the function to handles 

% Control algorithm
ctrlAlgorithm = @PIDControl;

% Simulation model
simModel = @MVPmodel;

% Observed variables
observationModel = @CHMsensor;

% Simulation method/function
simMethod = @ExplicitEuler;

% Objective function
objectiveFunction = @asymmetricQuadraticPenaltyFunction;

% Shooting function 
shootingmodel = @singleShootingObjective;

% Output model
outputModel = @CGMsensor;

%% Inizialization

% Scaling factor for the objective function (purely for numerical reasons)
scalingFactor = 1e-2;

% Index of the insulin bolus in the vector of manipulated inputs
idxbo = 2;

% Initial guess of the optimal insulin bolus
ubo0 = 0; % [mU/min]

N=15000/Ts;
M=100/Ts;

% Number of time steps in each control/sampling interval
NK=10;

% Meal and meal bolus after 1 hour
tMeal           = 1*h2min;          % [min]
idxMeal         = tMeal  /Ts + 1;   % [#]
D(1, idxMeal)   = 90     /Ts;       % [g CHO/min]

%% Simulation

% Optimal bolus 
ubo=computeOptimalBolus(N,M,idxbo,scalingFactor,x0,tspan,U,D,p,simModel,simMethod,objectiveFunction,NK);

% Meal and meal bolus  
U(idxbo, idxMeal) = ubo; 

% Simulate 
[T, X] = OpenLoopSimulation(x0, tspan, U, D, p, simModel, simMethod,NK);

% Blood glucose concentration 
G = CGMsensor(X, p); % [mg/dL]

%% Visualization
% Create figure with absolute size for reproducibility
figure;

% Plot blood glucose concentration
subplot(411);
plot(T*min2h, G);
xlim([t0, tf]*min2h);
ylabel({'Blood glucose concentration', '[mg/dL]'});

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*D(1, :), 'MarkerSize', 0.1);
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, U(1, [1:end, end]));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*U(2, :), 'MarkerSize', 1);
xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');





















