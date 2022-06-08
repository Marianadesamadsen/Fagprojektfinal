%% Simulation test of the GRID algorithm 30 days, 3 meals and two snacks while testing different values of Gmin and noise intensities.

% Simulating 3 meals on 30 days, and detecting the meals
% by using the GRID algorithm based on different values of Gmin and noise
% intensities.

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

%% Meal and meal bolus at 7, 12, 18 hours and snacks at 10, 15 hours

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

%% Making meal sizes and snack size. No bolus since goal is to detect all meals 

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

%% Simulating the control states based on x0 using different values of noise intensities. 

intensity_range = (1:10);                                           % Noise intensity range

% Initialising
T               = zeros(1,length(tspan),length(intensity_range));   % Three dimension. The third being for the different intensities.
X               = zeros(7,length(tspan),length(intensity_range));   % Three dimension. The third being for the different intensities.
G               = zeros(1,length(tspan),length(intensity_range));   % Three dimension. The third being for the different intensities.


% Looping over alle the values of intensities
for i = 1:length(intensity_range)
    
    [T(:,:,i), X(:,:,i)] = OpenLoopSimulation_withnoise(x0, tspan, U, D, p, @MVPmodel, @EulerM, Nk,intensity_range(i),5);
    
    % Extracting the blood glucose concentration 
    G(:,:,i) = CGMsensor_withnoise(X(:,:,i),p); % [mg/dL] 

end

%% Loading in the optimal values for Gmin based on the test from the deterministic simulation. 

Gmin_values       = load('Gminoptimal1patient.mat');
Gmin_combinations = Gmin_values.Gmin_optimal;           % All the optimal combination of Gmin values

%% GRID algorithm test of different intensities and Gmin values.

% Initialising
delta_G                      = 15;                                              % From article
t_vec                        = [5,10,15];                                       % The respective sampling times
tau                          = 6;                                               % From the article
intensity_number             = length(intensity_range);                         % Number of possibles intensity values
Gmin_number_combinations     = length(Gmin_combinations);                       % Number of combinations 
number_detectedmeals         = zeros(1,Gmin_number_combinations);               % Number of detected meals for each combination of Gmin

truepositive                 = zeros(intensity_number,Gmin_number_combinations);% True positives for each Gmin and intensity value
falsepositive                = zeros(intensity_number,Gmin_number_combinations);% False positives for each Gmin and intensity value
falsenegative                = zeros(intensity_number,Gmin_number_combinations);% False negative for each Gmin and intensity value
truenegative                 = zeros(intensity_number,Gmin_number_combinations);% True negative for each Gmin and intensity value

stride = 90/Ts; % How long it can possibly take to detect meal from the time the meal was given.

% Looping over possible intensities 
for j = 1 : intensity_number
    
    % Looping over all possible Gmin values for each possible intensity
    % value
    for i = 1 : Gmin_number_combinations
        
        % Detected meal for current Gmin value at current intensity
        D_detected_temp = GRIDalgorithm_mealdetection(G(:,:,j),Gmin_combinations(i,:),tau,delta_G,t_vec,Ts);
        
        % The total number of detected meals at current intensity and Gmin
        % value
        number_detectedmeals(j,i) = sum(D_detected_temp);
        
        % Computing all the true detected meals, false detected meals, and
        % not detected meals and true not detected meals. 
       [truenegative(j,i),truepositive(j,i),falsepositive(j,i),falsenegative(j,i)] = detectionrates(stride,D,D_detected_temp,Ts);

    end
    
end

% We can see we really fast get a lot of falsepositives when increasin the
% intensity

% We can see that we get more and more falsenegative slower than
% falsepositives when increasing the intensity


%% Calculating percentage error 

actualvalue = 90;

percentage = zeros(size(number_detectedmeals));

for j = 1 : intensity_number
    for i = 1 : Gmin_number_combinations
        percentage(j,i) = ( number_detectedmeals(j,i) - actualvalue ) / actualvalue;
    end
end

%% Calculating false positive and false negative rates 

for j = 1 : intensity_number

    falsepositive_rate(j,:) = falsepositive(j,:) ./ (falsepositive(j,:) + truenegative(j,:));
    truepositive_rate(j,:) = truepositive(j,:) ./ (truepositive(j,:) + falsenegative(j,:));

end

%% Visualize the ROC curve

figure 

subplot(2,5,1)
plot(falsepositive_rate(3,:),truepositive_rate(3,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 1')

subplot(2,5,2)
plot(falsepositive_rate(2,:),truepositive_rate(2,:),'*')
xlim([0 1])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 2')

subplot(2,5,3)
plot(falsepositive_rate(3,:),truepositive_rate(3,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 3')

subplot(2,5,4)
plot(falsepositive_rate(4,:),truepositive_rate(4,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 4')

subplot(2,5,5)
plot(falsepositive_rate(5,:),truepositive_rate(5,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 5')

subplot(2,5,6)
plot(falsepositive_rate(6,:),truepositive_rate(6,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 6')

subplot(2,5,7)
plot(falsepositive_rate(7,:),truepositive_rate(7,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 7')

subplot(2,5,8)
plot(falsepositive_rate(8,:),truepositive_rate(8,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 8')

subplot(2,5,9)
plot(falsepositive_rate(9,:),truepositive_rate(9,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 9')

subplot(2,5,10)
plot(falsepositive_rate(10,:),truepositive_rate(10,:),'*')
%xlim([0 0.00001])
ylim([0 1])
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC curve for intensity 10')










