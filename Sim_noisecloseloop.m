%% Sim closed loop
% Perform a closed-loop simulation with a PID controller with a 3 meals and 2 snacks over 30 days for the Medtronic
%%

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
ctrlAlgorithm = @PIDControl;

% Simulation model
simModel = @MVPmodel;

% Observed variables
observationModel = @CHMsensor;

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

% Index meals
idxMeal1         = tMeal1  /Ts + 1;   
idxMeal2         = tMeal2  /Ts + 1;   
idxMeal3         = tMeal3  /Ts + 1;   

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
        
end 

%% Simulate
% Closed-loop simulation
[T, X, Y, U] = ClosedLoopSimulation(tspan,x0,D,p, ... 
    ctrlAlgorithm, simMethod, simModel, observationModel, ctrlPar,ctrlState,Nk);

% Blood glucose concentration
Gsc = Y; % [mg/dL]

%% Detecting meals using GRID algorithm

% Inisializing 
filt_prev      = zeros(length(Gsc),2); % The vector of previous filtered values
Gfm_vec        = zeros(length(Gsc),2); % The vector of previous derivatives
G_vec          = [Gsc(1),Gsc(1),Gsc(1)];   % Inserting the start previous glucose measurements as the same value
delta_G        = 15;                 % From article
t_vec          = [5,10,15];          % The respective sampling times
filt_prev(1,:) = [Gsc(1),Gsc(1)];        % Inserting the previous filtered value as the not filtered values
tau            = 6;                  % From the article
flag           = 0;                  % No detected meal to begin with
Gmin           = [90 0.5 0.5];       % For meal under 50 considered

% Other tries
%Gmin = [ 130 1.5 1.6 ]; % Their meals
%Gmin = [ 110 1 1.5 ]; % For no meal under 50 


% Computing the first two detections
[ Gfm_vec(2,:) , filt_prev(2,:) , flag, zero_one(2) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(1,:) , Gmin, Gfm_vec(1,:) , t_vec ,flag );
                                
                                % Updating G_vec
                                G_vec=[Gsc(1),Gsc(1),Gsc(2)];
                                
[ Gfm_vec(3,:) , filt_prev(3,:) , flag, zero_one(3) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(2,:) , Gmin, Gfm_vec(2,:) , t_vec, flag );
                                
                                % Updating G_vec
                                G_vec=[Gsc(1),Gsc(2),Gsc(3)];
                                
% Computing the last detections
for i = 3 : length(Gsc)-1
    
[ Gfm_vec(i+1,:) , filt_prev(i+1,:) , flag, zero_one(i) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(i,:) , Gmin, Gfm_vec(i,:) , t_vec , flag );
                                
                                % Updating G_vec
                                G_vec=[Gsc(i-1),Gsc(i),Gsc(i+1)];
                                
end

% The total amount of detected meals
sum(zero_one)

%% Visualize
% Create figure with absolute size for reproducibility
figure;

% Plot blood glucose concentration
subplot(411);
plot(T*min2h, Gsc);
xlim([t0, tf]*min2h);
ylim([0 600])
ylabel({'CGM measurements', '[mg/dL]'});
hold on 
plot(tspan(1:end-1)*min2h,zero_one*200,'r.');

% Plot meal carbohydrate
subplot(412);
stem(tspan(1:end-1)*min2h, Ts*D(1, :), 'MarkerSize', 0.1);
xlim([t0, tf]*min2h);
ylim([-5 200])
ylabel({'Meal carbohydrates', '[g CHO]'});
hold on 
plot(tspan(1:end-1)*min2h,zero_one*100,'r.');

% Plot basal insulin flow rate
subplot(413);
stairs(tspan*min2h, U(1, [1:end, end]));
xlim([t0, tf]*min2h);
ylim([-5 100])
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan(1:end-1)*min2h, Ts*mU2U*U(2, :), 'MarkerSize', 1);
xlim([t0, tf]*min2h);
ylim([-1 1])
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');









