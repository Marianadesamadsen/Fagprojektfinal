%% Simulation test of the GRID algorithm 30 days, 3 meals with noise

% Simulating 3 meals on 30 days, and detecting the meals
% by using the GRID algorithm

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

%% Manipulated inputs
U = repmat(us, 1, N); % The same bolus and base rate for all

%% Disturbance variables
D = zeros(1, N); % No meal assumed

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

bolus = 0;
meal  = randi([50,150],1,90);

%% Inserting the meal sizes at the different hours/indicies

% Lopping over 30 days (one month)
for i = 0:29
    
    % Inserting the different meal sizes at the indcies 
        D(1, (idxMeal1+24*h2min/Ts*i))   = meal(1+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal1+24*h2min/Ts*i))   = bolus*U2mU/Ts;  
        D(1, (idxMeal2+24*h2min/Ts*i))   = meal(2+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal2+24*h2min/Ts*i))   = bolus*U2mU/Ts;  
        D(1, (idxMeal3+24*h2min/Ts*i))   = meal(3+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal3+24*h2min/Ts*i))   = bolus*U2mU/Ts;  
        
end 

%% Simulating the control states based on x0, the steady state.

[T, X] = OpenLoopSimulation(x0, tspan, U, D, p, @MVPmodel, @ExplicitEuler, Nk);

%% Blood glucose concentration 

G = CGMsensor(X, p); % [mg/dL] 

%% Detecting meals using GRID algorithm

% Inisializing 
filt_prev      = zeros(length(G),2); % The vector of previous filtered values
Gfm_vec        = zeros(length(G),2); % The vector of previous derivatives
G_vec          = [G(1),G(1),G(1)];   % Inserting the start previous glucose measurements as the same value
delta_G        = 15;                 % From article
t_vec          = [5,10,15];          % The respective sampling times
filt_prev(1,:) = [G(1),G(1)];        % Inserting the previous filtered value as the not filtered values
tau            = 6;                  % From the article
flag           = 0;                  % No detected meal to begin with
Gmin           = [120 0.8 0.65];     % Gmin accepts intensity up to 6

% Gmin = [90 0.5 0.5] % For no meal under 50 considered
% Other tries
% Gmin = [ 130 1.5 1.6 ]; % Their meals
%Gmin = [ 110 1 1.5 ]; % For no meal under 50 


% Computing the first two detections
[ Gfm_vec(2,:) , filt_prev(2,:) , flag, zero_one(2) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(1,:) , Gmin, Gfm_vec(1,:) , t_vec ,flag );
                                
                                % Updating G_vec
                                G_vec=[G(1),G(1),G(2)];
                                
[ Gfm_vec(3,:) , filt_prev(3,:) , flag, zero_one(3) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(2,:) , Gmin, Gfm_vec(2,:) , t_vec, flag );
                                
                                % Updating G_vec
                                G_vec=[G(1),G(2),G(3)];
                                
% Computing the last detections
for i = 3 : length(G)-1
    
[ Gfm_vec(i+1,:) , filt_prev(i+1,:) , flag, zero_one(i) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(i,:) , Gmin, Gfm_vec(i,:) , t_vec , flag );
                                
                                % Updating G_vec
                                G_vec=[G(i-1),G(i),G(i+1)];
                                
end

% The total amount of detected meals
sum(zero_one)

%% Visualize 

% Create figure with absolute size for reproducibility
figure;

% Plot blood glucose concentration and the detected meals as points
subplot(511);
plot(T*min2h, G);
xlim([t0, tf]*min2h);
ylabel({'Blood glucose concentration', '[mg/dL]'});
hold on 
plot(tspan(1:end-1)*min2h,zero_one*200,'r.');

% Plot meal carbohydrate and the detected meals as points
subplot(512);
stem(tspan(1:end-1)*min2h, Ts*D(1, :), 'MarkerSize', 0.1);
xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});
hold on 
plot(tspan(1:end-1)*min2h,zero_one*100,'r.');

% Plot basal insulin flow rate
subplot(513);
stairs(tspan*min2h, U(1, [1:end, end]));
xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(514);
stem(tspan(1:end-1)*min2h, Ts*mU2U*U(2, :), 'MarkerSize', 1);
xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[U]'}); 
xlabel('Time [h]');

% Plot detected Meals
subplot(515);
plot(tspan(1:end-1)*min2h,zero_one,'b-');
xlim([t0, tf]*min2h);
ylabel({'detected meal'}); 
xlabel('Time [h]'); 
 

