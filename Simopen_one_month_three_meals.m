%% Simulation openloop

% Simulate a single meal response in the MVP model
% by using x0 based on the steadystate vector

clear all 
clc 
close all 

%% Formatting the plots 

fs = 11; % Font size 
lw = 3; % Line width 
set(groot, 'DefaultAxesFontSize',   fs); % Set default font size 

% Set default line widths 
set(groot, 'DefaultLineLineWidth',  lw);
set(groot, 'DefaultStairLineWidth', lw);
set(groot, 'DefaultStemLineWidth',  lw);

% reset(groot);

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
tf = 720*h2min; % min - end time
Ts = 5;        % min - sampling time

%% Number of contral intervals

N = (tf - t0)/5; % [#]

%% Number of time steps in each control/sampling interval
Nk = 10;

%% Time span
tspan = 5*(0:N);

%% Initial condition
x0 = xs;

%% Manipulated inputs
U = repmat(us, 1, N); % The same bolus and base rate for all

%% Disturbance variables
D = zeros(1, N); % No meal assumed

%% Meal and meal bolus at hours 7,12,18
tMeal1           = 7*h2min;          % [min]
tMeal2           = 12*h2min;
tMeal3           = 18*h2min;
idxMeal1         = tMeal1  /Ts + 1;   % [#]
idxMeal2         = tMeal2  /Ts + 1;   % [#]
idxMeal3         = tMeal3  /Ts + 1;   % [#]

%% Making meal sizes with respectivly bolus sizes

meals=[50,70,10,120,40,80,110,90,60];
bolus=[6,8,2,12,5,9,12,10,7];

%% Inserting the meal sizes at the different hours/index

% Lopping over 30 days (one month)

for i = 0:29
    
    % Making new random number for the index of different meal sizes.
    k=randi(length(meals)-2);
    
    % Inserting the different meal sizes at the indcies 
        D(1, (idxMeal1+24*h2min/Ts*i))   = meals(k)     /Ts;       % [g CHO/min]
        U(2, (idxMeal1+24*h2min/Ts*i))   = bolus(k)*U2mU/Ts;  
        D(1, (idxMeal2+24*h2min/Ts*i))   = meals(k+1)     /Ts;       % [g CHO/min]
        U(2, (idxMeal2+24*h2min/Ts*i))   = bolus(k+1)*U2mU/Ts;  
        D(1, (idxMeal3+24*h2min/Ts*i))   = meals(k+2)     /Ts;       % [g CHO/min]
        U(2, (idxMeal3+24*h2min/Ts*i))   = bolus(k+2)*U2mU/Ts;  
end

%% Simulating the control states

[T, X] = OpenLoopSimulation(x0, tspan, U, D, p, @MVPmodel, @ExplicitEuler, Nk);

%% Blood glucose concentration 

G = CGMsensor(X, p); % [mg/dL] 

%% Visualize 

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
 

