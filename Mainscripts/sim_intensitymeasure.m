%% Closed loop simulation 
% This script is used to find the optimal intensity of the noise 

%%

clear all 
clc 
close all 

%% Loading all folders
fprintf('Loading diabetes library .. ');

% Add real thermodynamics functions
addpath(genpath(fullfile(pwd, '../Other')));

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
reset(groot)

%%  Conversion factors 

h2min = 60;      % Convert from h   to min 
min2h = 1/h2min; % Convert from min to h 
U2mU  = 1e3;     % Convert from U   to mU 
mU2U  = 1/U2mU;  % Convert from mU  to U 
min2sec = h2min; % Convert from min to sec
sec2min = 1/h2min;% Convert from sec to min

%% Inizializing parameters

p = [49,47,20.1,0.0106,0.0081,0.0022,1.33,253,47,5]; % Parameters at steady state
Gs = 108; % [mg/dL]: Steady state blood glucose concentration
ts = [];

%% Computing steadty state

[xs, us, flag] = computeSteadyStateMVPModel(ts, p, Gs);

% If fsolve did not converge, throw an error
if(flag ~= 1), error ('fsolve did not converge!'); end

%% Intital and final time for the simulation over 30 days

t0 =  0;        % min - start time
tf = 720*h2min; % min - end time
Ts = 5;         % min - step size 

%% Number of contral intervals

N = (tf - t0)/5; % [#]

%% Number of time steps in each control/sampling interval
Nk = 10;

%% Time span
tspan = 5*(0:N);

%% Initial condition
x0 = xs;

%% Inizialising the function to handles 
% Control algorithm
ctrlAlgorithm = @PIDControl2;

% Simulation model
simModel = @MVPmodel;

% Observed variables
observationModel = @CGMsensor_withnoise;

% Simulation method/function
simMethod = @EulerM;

%% Controller parameters and state
ctrlPar = [
      5.0;    % [min]     Sampling time
      0.05;   %           Proportional gain
      0.0005; %           Integral gain
      0.2500; %           Derivative gain
    108.0;    % [mg/dL]   Target blood glucose concentration
    NaN
    50.0;
    25.0;
    ];     % [mU/min]  Nominal basal rate (overwritten below)

ctrlState = [
      0.0;  %          Initial value of integral
    108.0]; % [mg/dL] Last measurements of glucose (dummy value)

%% Updating the nominal basal rate at steady state 

ctrlPar(6) = us(1);

%% Disturbance variables
D = zeros(1, N);

%% Meal and meal bolus at 7, 12, 18 hours

% Time meals
tMeal1           = 7*h2min;         
tMeal2           = 12*h2min;
tMeal3           = 18*h2min;
tSnack1          = 15*h2min;
tSnack2          = 10*h2min; 

% Index meals
idxMeal1         = tMeal1  /Ts + 1;   
idxMeal2         = tMeal2  /Ts + 1;   
idxMeal3         = tMeal3  /Ts + 1;
idxSnack1        = tSnack1 /Ts + 1;   
idxSnack2        = tSnack2 /Ts + 1;

%% Making meal sizes 

meal  = randi([50,150],1,90);
snack = 20;

%% Inserting the meal sizes at the different hours/indicies

% Lopping over 30 days (one month)
for i = 0:29
    
    % Inserting the different meal sizes at the indcies 
        D(1, (idxMeal1+24*h2min/Ts*i))   = meal(1+3*i)     /Ts;       % [g CHO/min]
        D(1, (idxMeal2+24*h2min/Ts*i))   = meal(2+3*i)     /Ts;       % [g CHO/min]
        D(1, (idxMeal3+24*h2min/Ts*i))   = meal(3+3*i)     /Ts;       % [g CHO/min]
       
    % Inserting the different meal sizes at the indcies 
        D(1, (idxSnack1+24*h2min/Ts*i))   = snack        /Ts;       % [g CHO/min] 
        D(1, (idxSnack2+24*h2min/Ts*i))   = snack        /Ts;       % [g CHO/min]
        
end 

%% 

% Insulin vector 
U = zeros(2,length(D(1,:)));

% Calculating the insulin amount for each meal 
for i = 1 : length(U)
    if D(1,i) > 50/Ts %because snack
        U(2,i) = (D(1,i)*Ts/10)*U2mU/Ts;
    end
end

% Removing bolus insulin from some of the meal indices. 
% Every 5th day the meal at 7 hour is missed.
for i = 1 : 5 : 29 
    U(2,idxMeal1+24*h2min/Ts*i) = 0;
end

% Decreasing the amount of insulin for some of the meal indices. 
% Every 5th day starting at day 3 the bolus insulin is 0.5 too low. 
for i = 3 : 5 : 29
    U(2,idxMeal2+24*h2min/Ts*i) = U(2,idxMeal2+24*h2min/Ts*i) * 0.5;
end


%% Simulate and visualize

% Create figure with absolute size for reproducibility
figure;

i = 1;

for k = 0:5:20

intensity = k;

% Closed-loop simulation
[T, X, Y, U] = ClosedLoopSimulation_withnoise2(tspan,x0,D,U,p, ... 
    ctrlAlgorithm, simMethod, simModel, observationModel, ctrlPar,ctrlState,Nk,intensity);

% Converting data for axis
T2=datetime(T*min2sec,'ConvertFrom','posixtime');
tspan2=datetime(tspan*min2sec,'ConvertFrom','posixtime');

% Blood glucose concentration
Gsc = Y; % [mg/dL]

% Plot blood glucose concentration
subplot(length(0:5:20),1,i);
plot(T2(1:720), Gsc(1:720)); %visualizing three day 
%xlim([t0, tf]*min2h);
ylim([0 400])
ylabel({'CGM measurements', '[mg/dL]'});
title(sprintf('Intensity %d',k))

i = i+1;

end
