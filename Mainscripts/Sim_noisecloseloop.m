%% Sim closed loop
% Perform a closed-loop simulation with a PID controller with a 3 meals
% and 2 snacks over 30 days for the Medtronic with stochastic noise
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

p = [49,47,20.1,0.0106,0.0081,0.0022,1.33,253,47,5]; % Parameters at steady state
Gs = 108; % [mg/dL]: Steady state blood glucose concentration
ts = [];

%% Computing steady state

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

%% Removing bolus insulin every 5th day at 7 AM and decreasing bolus insulin every 5th day starting from day 3 at 7 AM

% Insulin vector
U = zeros(2,length(D(1,:)));
idx_missed_temp = zeros(1,26);
idx_less_temp = zeros(1,28);

% Calculating the insulin amount for each meal
for i = 1 : length(U)
    if D(1,i) > 50/Ts %because snack
        U(2,i) = (D(1,i)*Ts/10)*U2mU/Ts; % ICR
    end
end

% Removing bolus insulin from some of the meal indices.
% Every 5th day the meal at 7 hour is missed.
for i = 1 : 3 : 29
    U(2,idxMeal1+24*h2min/Ts*i) = 0;
    idx_missed_temp(i) = idxMeal1+24*h2min/Ts*i ;
end

idx_missed = nonzeros(idx_missed_temp)';

% Decreasing the amount of insulin for some of the meal indices.
% Every 5th day starting at day 3 the bolus insulin is 0.5 too low.
for i = 2 : 3 : 29
    U(2,idxMeal2+24*h2min/Ts*i) = U(2,idxMeal2+24*h2min/Ts*i) * 0.2;
    idx_less_temp(i) = idxMeal2+24*h2min/Ts*i;
end

idx_less = nonzeros(idx_less_temp)';


%% Simulate

intensity = 10;

% Closed-loop simulation
[T, X, Y, U] = ClosedLoopSimulation_withnoise2(tspan,x0,D,U,p, ...
    ctrlAlgorithm, simMethod, simModel, observationModel, ctrlPar,ctrlState,Nk,intensity);

% Blood glucose concentration
Gsc = Y; % [mg/dL]

%% Compting the different combinations of Gmin value based on their different rates

% Range for gmin try outs
gmin1range = (130:5:140);
gmin2range = (1.5:0.5:2.5);
gmin3range = (1.3:0.5:2.3);

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
number_detectedmeals    = zeros(1,number_combinations);
truepositive            = zeros(1,number_combinations);
falsepositive           = zeros(1,number_combinations);
falsenegative           = zeros(1,number_combinations);
truenegative            = zeros(1,number_combinations);
D_detected              = zeros(1,N,number_combinations);

stride = 90/Ts; % How long it can possibly take to detect meal from the time the meal was given.

% Looping over all the different combinations of Gmin values
for i = 1 : number_combinations(1)

% Detecting meals
D_detected(:,:,i) = GRIDalgorithm_mealdetection(Gsc,Gmin_combinations(i,:),tau,delta_G,t_vec,Ts);

% Total number of detected meals for the current Gmin values.
number_detectedmeals(i) = sum(D_detected(:,:,i));

[truepositive(i), falsepositive(i), falsenegative(i), truenegative(i)] = detectionrates2(stride,D,D_detected(:,:,i),Ts,idx_missed, idx_less,U);

end

Detec = number_detectedmeals';
TP = truepositive';
FP = falsepositive';
FN = falsenegative';
TN = truenegative';

table(Detec,TP,FP,FN,TN)

%% Computing the rates to find the optimal Gmin values

rateFP = FP/30; % The mean of false positives pr day
rateTP = TP/20; % The percent of how many true positives out of the total
idxFPbest = 1;

for i = 1 : number_combinations  
    
    % We want to have at max 0.5 false positives pr day
    if rateFP(i) <= 0.5
        idxFPbest = i;
    end
    
    % We want to detect at least 70% of the meals we want to detect
    % We want to have at max 0.5 false positives pr day
    if rateFP(idxFPbest) <= 0.5 && rateTP(idxFPbest) >= 0.7
        idx_tempoptimal(i) = i;
    end
      
end

idx_optimalfinal = nonzeros(idx_tempoptimal');

rFP_optimal = rateFP(idx_optimalfinal);
rTP_optimal = rateTP(idx_optimalfinal);
rFP_total = rateFP;
rTP_total = rateTP;

table(rFP_optimal,rTP_optimal,idx_optimalfinal)
table(rFP_total,rTP_total)

%% Vectors with length of time of when bolus is missed and lessened

% Converting data
T2=datetime(T*min2sec,'ConvertFrom','posixtime');
tspan2=datetime(tspan*min2sec,'ConvertFrom','posixtime');

missed_vector = zeros(1,length(T2)-1);
for i = 1:length(idx_missed)
    k = idx_missed(i);
    missed_vector(k) = 1;
end

less_vector = zeros(1,length(T2)-1);
for i = 1:length(idx_less)
    k = idx_less(i);
    less_vector(k) = 1;
end


%% Visualize
% Create figure with absolute size for reproducibility
figure;
% Plot blood glucose concentration
subplot(411);
plot(T2, Gsc);
%xlim([t0, tf]*min2h);
%ylim([0 300])
ylabel({'CGM measurements', '[mg/dL]'});
title('Blood glucose concentration over time')
hold on
plot(tspan2(1:end-1),D_detected(:,:,end-7)*200,'r.');
hold on
plot(tspan2(1:end-1),missed_vector*150,'b *');
hold on
plot(tspan2(1:end-1),less_vector*150,'g *');

% Plot meal carbohydrate
subplot(412);
stem(tspan2(1:end-1), Ts*D(1, :), 'MarkerSize', 0.1);
%xlim([t0, tf]*min2h);
ylim([-5 200])
ylabel({'Meal carbohydrates', '[g CHO]'});
%hold on
title('Meals and meal sizes')
%plot(tspan2(1:end-1),D_detected*100,'r.');
hold on
plot(tspan2(1:end-1),missed_vector*50,'b *');
hold on
plot(tspan2(1:end-1),less_vector*50,'g *');

% Plot basal insulin flow rate
subplot(413);
stairs(tspan2, U(1, [1:end, end]));
%xlim([t0, tf]*min2h);
ylim([-5 100])
ylabel({'Basal insulin', '[mU/min]'});
title('Basal insulin flow rate')

% Plot bolus insulin
subplot(414);
stem(tspan2(1:end-1),Ts*mU2U*U(2, :), 'MarkerSize', 1);
%xlim([t0, tf]*min2h);
%ylim([0 5])
ylabel({'Bolus insulin', '[U]'});
xlabel('Time [h]');
title('Bolus insulin')
