%% Simulation test of the GRID algorithm on 100 patients

% Simulating 3 meals and 2 snacks on 30 days on 100 patients at a time
% using the GRID algorithm

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

numpatients = 100;         % number of patients
pf = pmatrix(numpatients); % computing the p vectors for all patients
Gs = 108;                  % Steady state blood glucose concentration
ts = [];

%% Computing steadty state

% Inisializing
xs=zeros(numpatients,:);
us=zeros(numpatients,:);

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
U=zeros(:,:,numpatients);

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

% meals=[50,70,10,120,40,80,110,90,60];
% meals=[50,70,120,80,110,90,60];
% bolus=[6,8,2,12,5,9,12,10,7];
% bolus=zeros(1,length(meals)); % Bolus is zero because we will detect meals

bolus = 0;
meal  = randi([50,150],1,90);
snack = 20;

%% Inserting the meal sizes at the different hours/index

% Lopping over 30 days (one month)

for p=1:numpatients

for i = 0:29
    
    % Inserting the different meal sizes at the indcies 
        D(1, (idxMeal1+24*h2min/Ts*i),p)   = meal(1+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal1+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        D(1, (idxMeal2+24*h2min/Ts*i),p)   = meal(2+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal2+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        D(1, (idxMeal3+24*h2min/Ts*i),p)   = meal(3+3*i)     /Ts;       % [g CHO/min]
        U(2, (idxMeal3+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        
    % Inserting the different meal sizes at the indcies 
        D(1, (idxSnack1+24*h2min/Ts*i),p)   = snack     /Ts;       % [g CHO/min]
        U(2, (idxSnack1+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        D(1, (idxSnack2+24*h2min/Ts*i),p)   = snack    /Ts;       % [g CHO/min]
        U(2, (idxSnack2+24*h2min/Ts*i),p)   = bolus*U2mU/Ts;  
        
end

end

%% Simulating the control states

for p=1:numpatients

[T(:,:,p), X(:,:,p)] = OpenLoopSimulation(x0(p,:)', tspan, U(:,:,p), D(:,:,p), pf(:,p), @MVPmodel, @ExplicitEuler, Nk);

end

%% Blood glucose concentration 

for p=1:numpatients

G(:,:,p) = CGMsensor(X(:,:,p), pf(:,p)); % [mg/dL] 

end

%% GRID

prev_vec = zeros(length(G(:,:,1)),2,numpatients);
Gf_vec = zeros(length(G),2,numpatients);
detectedmeals = zeros(1,numpatients);

for p = 1:numpatients 

G_grid = [G(:,1,p),G(:,1,p),G(:,1,p)];
delta_G = 15;
tspan2 = 5;
t_vec = [5,10,15];
prev_vec(1,:,p) = [G(:,1,p),G(:,1,p)];
% Gmin = [ 130 1.5 1.6 ]; % Their meals
% Gmin = [ 110 1 1.5 ]; % For no meal under 50 
Gmin = [90 0.5 0.5]; % For meal under 50 considered

tau = 6; 
flag = 0; 

[ Gf_vec(2,:,p) , prev_vec(2,:,p) , flag, zero_one(2) ] = GRID_func( delta_G , G_grid , tau, tspan2 , ...
                                    prev_vec(1,:,p) , Gmin, Gf_vec(1,:,p) , t_vec ,flag );
                                
                                G_grid=[G(:,1,p),G(:,1,p),G(:,2,p)];
                                
[ Gf_vec(3,:,p) , prev_vec(3,:,p) , flag, zero_one(3) ] = GRID_func( delta_G , G_grid , tau, tspan2 , ...
                                    prev_vec(2,:,p) , Gmin, Gf_vec(2,:,p) , t_vec, flag );
                                
                                G_grid=[G(:,1,p),G(:,2,p),G(:,3,p)];

for i = 3 : length(G)-1
    
[ Gf_vec(i+1,:,p) , prev_vec(i+1,:,p) , flag, zero_one(i,p) ] = GRID_func( delta_G , G_grid , tau, tspan2 , ...
                                    prev_vec(i,:,p) , Gmin, Gf_vec(i,:,p) , t_vec , flag );
                                
                                G_grid=[G(i-1),G(i),G(i+1)];
                                
end

detectedmeals(p)=sum(zero_one(:,p));

end

