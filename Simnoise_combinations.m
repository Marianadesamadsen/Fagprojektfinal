%% Simulation test of the GRID algorithm 30 days, 3 meals and two snacks with noise

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

%% Simulating the control states based on x0, the steady state over different types of intensities

intensity_range = (1:10);
T         = zeros(1,length(tspan),length(intensity_range));
X         = zeros(7,length(tspan),length(intensity_range));
G         = zeros(1,length(tspan),length(intensity_range));

for i = 1:length(intensity_range)
    
    [T(:,:,i), X(:,:,i)] = OpenLoopSimulation_withnoise(x0, tspan, U, D, p, @MVPmodel, @EulerM, Nk,intensity_range(i));
    
    % Blood glucose concentration 
    G(:,:,i) = CGMsensor_withnoise(X(:,:,i),p); % [mg/dL] 

end

%% Loading in the optimal values for Gmin

Gmin_values       = load('Gminoptimal1patient.mat');
Gmin_combinations = Gmin_values.Gmin_optimal;

%%

% Initialising
delta_G        = 15;                 % From article
t_vec          = [5,10,15];          % The respective sampling times
tau            = 6;                  % From the article

Gmin_number_combinations     = size(Gmin_combinations); 
number_detectedmeals         = zeros(1,Gmin_number_combinations(1));
truepositive                 = zeros(1,Gmin_number_combinations(1));
falsepositive                = zeros(1,Gmin_number_combinations(1));
falsenegative                = zeros(1,Gmin_number_combinations(1));
truenegative                 = zeros(1,Gmin_number_combinations(1));

stride = 90/Ts; % min


for j = 1 : length(intensity_range)
    
    for i = 1 : Gmin_number_combinations(1)

        D_detected = GRIDalgorithm_mealdetection(G(:,:,j),Gmin_combinations(i,:),tau,delta_G,t_vec,Ts);

        number_detectedmeals(j,i) = sum(D_detected);

        [truenegative(j,i),truepositive(j,i),falsepositive(j,i),falsenegative(j,i)] = detectionrates(stride,D,D_detected,Ts);

    end
    
end


%% Calculating false positive and false negative rates 

falsepositive_rate = zeros(length(intensity_range),Gmin_number_combinations(1));
truepositive_rate  = zeros(length(intensity_range),Gmin_number_combinations(1));

for j = 1 : length(intensity_range)

falsepositive_rate(j,:) = falsepositive(j,:) ./ (falsepositive(j,:) + truenegative(j,:));
truepositive_rate(j,:)  = truepositive(j,:) ./ (truepositive(j,:) + falsenegative(j,:));

end

%% Visualize 

figure 

plot(falsepositive_rate(1,:),truepositive_rate(1,:),'*')
xlim([0 1.2*10^(-4)])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve')

%% Finding the optimal index for intensity and gmin 

gmin_idx        = zeros(1,length(falsepositive_rate));
intensity_idx   = zeros(1,length(intensity_range));

for j = 1 : length(intensity_range)
    
    for i = 1 : length(falsepositive_rate)
   
        if falsepositive_rate(j,i) == 0 && truepositive_rate(j,i) == 1
            gmin_idx(i)      = i; 
            intensity_idx(j) = j;
        end 
        
    end
    
end

gmin_idx_best     = find(gmin_idx);
itensity_idx_best = find(gmin_idx);

%%

Gmin_optimal = zeros(length(gmin_idx_best),3);

for i = 1 : length(gmin_idx_best)
    
   k=gmin_idx_best(i);
   
   Gmin_optimal(i,:) = Gmin_combinations(k,:);
        
end











