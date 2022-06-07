%% Simulation test of the GRID algorithm on 100 patients

% Simulating 3 meals and 2 snacks on 30 days on 100 patients at a time
% using the GRID algorithm

%% 

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

%%  Conversion factors 

h2min = 60;      % Convert from h   to min 
min2h = 1/h2min; % Convert from min to h 
U2mU  = 1e3;     % Convert from U   to mU 
mU2U  = 1/U2mU;  % Convert from mU  to U
min2sec = h2min; % Convert from min to sec
sec2min = 1/h2min;% Convert from sec to min

%% Inizializing parameters

numpatients = 100;         % number of patients
pf = pmatrix(numpatients); % computing the p vectors for all patients
Gs = 108;                  % Steady state blood glucose concentration
ts = [];

%% Computing steadty state

% Inisializing
xs=zeros(numpatients,7);
us=zeros(numpatients,2);

% Looping over all patients
for i=1:numpatients
    
    [xs(i,:), us(i,:), flag] = computeSteadyStateMVPModel(ts, pf(:,i), Gs);

    % If fsolve did not converge, throw an error
    if(flag ~= 1), error ('fsolve did not converge!'); end

end

%% Intital and final time for the simulation over 30 days
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

% Inisializing
U=zeros(2,N,numpatients);

% Looping over all patients 
for i=1:numpatients
    
    U(:,:,i) = repmat(us(i,:)', 1, N); % The same bolus and base rate for all time samples

end

%% Disturbance variables
D = zeros(1, N,numpatients); % No meal assumed

%% Meals and snacks at 7,12,18 hours and 10,15 hours

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

%% Making meal sizes with respectivly bolus sizes

bolus = 0;
meal  = randi([50,150],1,90);
snack = 20;

%% Inserting the meal sizes at the different hours/indicies

% Looping over all patients
for p=1:numpatients

    % Looping over 30 days (one month)
    for i = 0:29
    
    % Inserting the different meal sizes at the indcies 
        D(1, (idxMeal1+24*h2min/Ts*i),p)   = meal(1+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal1+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        D(1, (idxMeal2+24*h2min/Ts*i),p)   = meal(2+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal2+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        D(1, (idxMeal3+24*h2min/Ts*i),p)   = meal(3+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal3+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        
    % Inserting the different meal sizes at the indcies 
        D(1, (idxSnack1+24*h2min/Ts*i),p)   = snack     /Ts;            % [g CHO/min]
        U(2, (idxSnack1+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        D(1, (idxSnack2+24*h2min/Ts*i),p)   = snack    /Ts;             % [g CHO/min]
        U(2, (idxSnack2+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        
     end

end

%% Simulating the control states for all patients

% Inisializing 
T=zeros(1,N+1,numpatients);
X=zeros(7,N+1,numpatients);

% Looping over all patients
for p=1:numpatients

[T(:,:,p), X(:,:,p)] = OpenLoopSimulation(x0(p,:)', tspan, U(:,:,p), D(:,:,p), pf(:,p), @MVPmodel, @ExplicitEuler, Nk);

end

%% Blood glucose concentration 

% Inisializing
G=zeros(1,N+1,numpatients);

% Looping over all patients
for p=1:numpatients

G(:,:,p) = CGMsensor(X(:,:,p), pf(:,p)); % [mg/dL] 

end

%% GRID

% Inisializing
number_detectedmeals  = zeros(1,numpatients);            
D_detected            = zeros(length(G)-1,numpatients);     

for p = 1:numpatients 

% Initializing
G_vec            = G(1,:,p);     % The current person
Gmin             = [90 0.5 0.5]; % Temperaly chosen
tau              = 6;            % From article   
delta_G          = 15;           % From article   
t_vec            = [5,10,15];    % The time is the same for each
  
% Computing the detected meals
D_detected(:,p) = GRIDalgorithm_mealdetection(G,Gmin,tau,delta_G,t_vec,Ts);

% The total amount of detected meals for each patient in vector
number_detectedmeals(p)=sum(D_detected(:,p));

end

