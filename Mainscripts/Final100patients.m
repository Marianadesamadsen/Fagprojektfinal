%% Closed loop simulation 100 patients
% Perform a closed-loop simulation with a PID controller with a 3 meals
% and 2 snacks over 30 days for the Medtronic with stochastic noise
% Detecting meals using the GRID algorithm with 125 combinations of the
% gmin values.
% Computes the optimal Gmin values using the evaluation function.

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

%% Initializing parameters

numpatients = 100;         % number of patients
pf = pmatrix(numpatients); % computing the p vectors for all patients
Gs = 108;                  % [mg/dL]: Steady state blood glucose concentration
ts = [];

%% Computing steadty state

% Initializing
xs=zeros(numpatients,7);
us=zeros(numpatients,2);

% Looping over all patients
for p=1:numpatients
    
    [xs(p,:), us(p,:), flag] = computeSteadyStateMVPModel(ts, pf(:,p), Gs);

    % If fsolve did not converge, throw an error
    if(flag ~= 1), error ('fsolve did not converge!'); end

end

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
D = zeros(N,numpatients);

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

% Looping over patients 
for p = 1 : numpatients
    
% Looping over 30 days (one month)
    for i = 0:29

    % Inserting the different meal sizes at the indcies
        D((idxMeal1+24*h2min/Ts*i),p)   = meal(1+3*i)     /Ts;       % [g CHO/min]
        D((idxMeal2+24*h2min/Ts*i),p)   = meal(2+3*i)     /Ts;       % [g CHO/min]
        D((idxMeal3+24*h2min/Ts*i),p)   = meal(3+3*i)     /Ts;       % [g CHO/min]

    % Inserting the different meal sizes at the indcies
        D((idxSnack1+24*h2min/Ts*i),p)   = snack        /Ts;       % [g CHO/min]
        D((idxSnack2+24*h2min/Ts*i),p)   = snack        /Ts;       % [g CHO/min]

    end
end

%% Removing bolus insulin every 5th day at 7 AM and decreasing bolus insulin every 5th day starting from day 3 at 7 AM

% Insulin vector
U = zeros(2,N,numpatients);

% Looping over patients
for p = 1 : numpatients

    % Calculating the insulin amount for each meal
    for i = 1 : length(U)
        if D(i,p) > 50/Ts %because snack
            U(2,i,p) = (D(i,p)*Ts/10)*U2mU/Ts; % ICR
        end
    end

    % Removing bolus insulin from some of the meal indices.
    % Every 5th day the meal at 7 hour is missed.
    for i = 1 : 3 : 29
        U(2,idxMeal1+24*h2min/Ts*i,p) = 0;
    end

    % Decreasing the amount of insulin for some of the meal indices.
    % Every 5th day starting at day 3 the bolus insulin is 0.5 too low.
    for i = 2 : 3 : 29
        U(2,idxMeal2+24*h2min/Ts*i,p) = U(2,idxMeal2+24*h2min/Ts*i,p) * 0.2;
    end
    
end

%% Simulate

intensity = 10;
% Inisializing 
T         = zeros(N+1,numpatients); 
X         = zeros(7,N+1,numpatients); 
Y         = zeros(N+1,numpatients);
Gsc_temp  = zeros(N+1,numpatients);

for p = 1 : numpatients

    % Closed-loop simulation
    [T(:,numpatients), X(:,:,p), Y(:,p), U(:,:,p)] = ClosedLoopSimulation_withnoise2(tspan,x0(p,:)',D(:,p)',U(:,:,p),pf(:,p), ...
        ctrlAlgorithm, simMethod, simModel, observationModel, ctrlPar,ctrlState,Nk,intensity);

    % Making sure to remove unusual patients
    if max(Y(:,p)) < 550 && min(Y(:,p)) > 20
    
        % Blood glucose concentration
        Gsc_temp(:,p) = Y(:,p); % [mg/dL]
        idxpatientstemp(p) = p;
    
    end
    
end

%% Number of new patients

idxpatientsfinal = nonzeros(idxpatientstemp');
numpatients_new = length(idxpatientsfinal);


%% Storing new glucose values to remove zeros 

% Initialising 
Gsc = zeros(N+1,numpatients_new);

Gsc_temp2 = nonzeros(Gsc_temp);

for p = 0:numpatients_new - 1 
    
    Gsc(:,p+1) = Gsc_temp2(1+p*(N+1):(1+p)*(N+1));
    
end

%% Compting the different combinations of Gmin value based on their different rates

% Range for gmin try outs
gmin1range = (120:5:140);
gmin2range = [0.6,1.1,1.6,2.1,2.6];
gmin3range = [0.5,1.0,1.5,2.0,2.5];

% Computing the combinations using meshgrid
[Gmin1,Gmin2,Gmin3] = meshgrid(gmin1range,gmin2range,gmin3range);

% All types of combinations
Gmin_combinations = [Gmin1(:),Gmin2(:),Gmin3(:)];

%% Detecting meals using GRID algorithm

% Inisializing
delta_G        = 15;                 % From article
t_vec          = [5,10,15];          % The respective sampling times
tau            = 6;                  % From the article

% Initialising
number_combinations     = length(Gmin_combinations);
number_detectedmeals    = zeros(number_combinations,numpatients_new);
truepositive            = zeros(number_combinations,numpatients_new);
falsepositive           = zeros(number_combinations,numpatients_new);
D_detected              = zeros(N,number_combinations,numpatients_new);

stride = 90/Ts; % How long it can possibly take to detect meal from the time the meal was given.

for p = 1 : numpatients_new 
    
    % Looping over all the different combinations of Gmin values
    for i = 1 : number_combinations(1)

    % Detecting meals
    D_detected(:,i,p) = GRIDalgorithm_mealdetection(Gsc(:,p)',Gmin_combinations(i,:),tau,delta_G,t_vec,Ts);

    % Total number of detected meals for the current Gmin values.
    number_detectedmeals(i,p) = sum(D_detected(:,i,p));

    [truepositive(i,p),falsepositive(i,p)] = evaluationfunction(stride,D(:,p)',D_detected(:,i,p),Ts,N);

    end

end

%% Computing the mean of all truenegatives, truepositives, falsepositives, falsenegatives for each combination

% Initializing 
meanFN = zeros(1,number_combinations);
meanTP = zeros(1,number_combinations);
meanDetec = zeros(1,number_combinations);

% Looping over all combinations
for i = 1 : number_combinations
    meanFP(i) = mean(falsepositive(i,:));
    meanTP(i) = mean(truepositive(i,:));
    meanDetec(i) = mean(number_detectedmeals(i,:));
end

%% Making a table 

MeanTP = meanTP';
MeanFP = meanFP';
MeanDetec = meanDetec';

%% Computing the rates to find the optimal Gmin values

rateFP = meanFP/30; % The mean of false positives pr day
rateTP = meanTP/90; % The percent of how many true positives out of the total
idx_tempoptimal = zeros(1,length(number_combinations));
for i = 1 : number_combinations 
    
    % We want to detect at least 70% of the meals we want to detect
    % We want to have at max 0.5 false positives pr day
    if rateFP(i) <= 1 && rateTP(i) >= 0.7
        idx_tempoptimal(i) = i;
    end
      
end

%%

idx_optimalfinal = nonzeros(idx_tempoptimal');

rFP_optimal = rateFP(idx_optimalfinal)';
rTP_optimal = rateTP(idx_optimalfinal)';
rFP_total = rateFP';
rTP_total = rateTP';

table(rFP_total,rTP_total)
table(rFP_optimal,rTP_optimal,idx_optimalfinal)

%%

close all

% Converting data
T2=datetime(T*min2sec,'ConvertFrom','posixtime');
tspan2=datetime(tspan*min2sec,'ConvertFrom','posixtime');

% Plot blood glucose concentration and the detected meals as points
subplot(411);
plot(T2, Y(:,32))
ylabel({'Blood glucose concentration', '[mg/dL]'});
title('Blood glucose concentration over time','FontSize', 25)

% Plot meal carbohydrate and the detected meals as points
subplot(412);
stem(tspan2(1:end-1), Ts*D(:,32)', 'MarkerSize', 0.1)
%xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});
title('Meals and meal sizes','FontSize', 25)

% Plot basal insulin flow rate
subplot(413);
stairs(tspan2, U(1, [1:end, end],32))
%xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});
title('Basal insulin flow rate','FontSize', 25)

% Plot bolus insulin
subplot(414);
stem(tspan2(1:end-1), Ts*mU2U*U(2, :,32), 'MarkerSize', 1)
%xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[U]'}); 
xlabel('Time [h]');
title('Bolus insulin with missing bolus and less bolus at some meals','FontSize', 25)



