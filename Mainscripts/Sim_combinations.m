%% Simulation test of the GRID algorithm 30 days, 3 meals and 2 snacks with different value combinations of Gmin

% Simulating 3 meals on 30 days, and detecting the meals
% by using the GRID algorithm with different value combinations of Gmin.

%%
clear all 
clc 

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

%%  Conversion factors 

h2min = 60;      % Convert from h   to min 
min2h = 1/h2min; % Convert from min to h 
U2mU  = 1e3;     % Convert from U   to mU 
mU2U  = 1/U2mU;  % Convert from mU  to U 
min2sec = h2min; % Convert from min to sec
sec2min = 1/h2min; %Convert from sec to min

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

%% Number of control intervals

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
tSnack1          = 15*h2min;
tSnack2          = 10*h2min;  

% Index meals
idxMeal1         = tMeal1  /Ts + 1;   
idxMeal2         = tMeal2  /Ts + 1;   
idxMeal3         = tMeal3  /Ts + 1;   
idxSnack1        = tSnack1 /Ts + 1;   
idxSnack2        = tSnack2 /Ts + 1;

%% Making meal sizes 

bolus = 0;
meal  = randi([50,150],1,90);
snack = 20; 

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
        
         % Inserting the different meal sizes at the indcies 
        D(1, (idxSnack1+24*h2min/Ts*i))   = snack     /Ts;            % [g CHO/min]
        U(2, (idxSnack1+24*h2min/Ts*i))   = bolus*U2mU/Ts;  
        D(1, (idxSnack2+24*h2min/Ts*i))   = snack    /Ts;             % [g CHO/min]
        U(2, (idxSnack2+24*h2min/Ts*i))   = bolus*U2mU/Ts;  
        
end 

%% Simulating the control states based on x0, the steady state.

[T, X] = OpenLoopSimulation(x0, tspan, U, D, p, @MVPmodel, @ExplicitEuler, Nk);

%% Blood glucose concentration 

G = CGMsensor(X, p); % [mg/dL] 

%% Compting the different combinations of Gmin value based on their different rates

% Range for gmin try outs 
gmin1range = (125:5:135);
gmin2range = (1.2:0.4:2.2);
gmin3range = (1.2:0.4:2.2);

% Computing the combinations using meshgrid
[Gmin1,Gmin2,Gmin3] = meshgrid(gmin1range,gmin2range,gmin3range);

% All types of combinations
Gmin_combinations = [Gmin1(:),Gmin2(:),Gmin3(:)];

%% Detecting meals using GRID algorithm and computing the number of FP,TP,TN,FN

% Inisializing 
delta_G        = 15;                 % From article
t_vec          = [5,10,15];          % The respective sampling times
tau            = 6;                  % From the article

% Initialising
number_combinations     = length(Gmin_combinations); 
number_detectedmeals    = zeros(1,number_combinations);
truepositive            = zeros(1,number_combinations);
falsepositive           = zeros(1,number_combinations);
falsenegative           = zeros(1,number_combinations);
truenegative            = zeros(1,number_combinations);

stride = 90/Ts; % How long it can possibly take to detect meal from the time the meal was given.

% Looping over all the different combinations of Gmin values
for i = 1 : number_combinations(1)

% Detecting meals
D_detected = GRIDalgorithm_mealdetection(G,Gmin_combinations(i,:),tau,delta_G,t_vec,Ts);

% Total number of detected meals for the current Gmin values.
number_detectedmeals(i) = sum(D_detected);

[truenegative(i),truepositive(i),falsepositive(i),falsenegative(i)] = detectionrates(stride,D,D_detected,Ts);

end

%% Calculating false positive and false negative rates 

falsepositive_rate = falsepositive ./ (falsepositive + truenegative);
truepositive_rate = truepositive ./ (truepositive + falsenegative);

%% Visualize the ROC curve

figure 

plot(falsepositive_rate,truepositive_rate,'*')
xlim([0 0.00001])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve')

%% Finding the indices for the optimal Gmin value, such that AUC = 1, meaning FPR=0 and TPR=1.

gmin_idx=zeros(1,length(falsepositive_rate));

for i = 1 : length(falsepositive_rate)
   
    if falsepositive_rate(i) == 0 && truepositive_rate(i) == 1
        gmin_idx(i) = i; 
    end 
    
end

% Indices for optimal gmin values
gmin_idx_optimal=find(gmin_idx);

%% Extracing the optimal Gmin values. 

Gmin_optimal = zeros(length(gmin_idx_optimal),3);

for i = 1 : length(gmin_idx_optimal)
    
   k=gmin_idx_optimal(i);
   
   Gmin_optimal(i,:) = Gmin_combinations(k,:);
        
end

% Saving the values. 

% It is in a comment since it is already done once, and
% will change for every time it is being runned. 

% save('Gminoptimal1patient','Gmin_optimal');









