%% Loading real clinical data from excel where data is stored in sheets
clear all
clc
%% Adding the path 
addpath(genpath(fullfile(pwd, '../Other')));

CGM_sheet = readtable('Control-IQ_Sample_Diasend','Sheet','CGM');
Insulin_carbs_sheet = readtable('Control-IQ_Sample_Diasend', 'Sheet','Insulin and Carbs');

% Store the data in temporary arrays and the time in cell
CGM = table2array(CGM_sheet(:,2));                  % Glucose concentration
uba_temp = table2array(Insulin_carbs_sheet(:,2));   % Basal insulin
ubo_temp = table2array(Insulin_carbs_sheet(:,4));   % Bolus insulin
meal_temp = table2array(Insulin_carbs_sheet(:,8));  % Amount of carbs

time_temp_CGM = table2cell(CGM_sheet(:,1));         % Time frame for glucose measurements
time_temp_insulin_carb = table2cell(Insulin_carbs_sheet(:,1)); % Time frame for insulin and meal measurements

%% Replacing the missing values NaN with zeros

% Replacing all NaN values with 0
uba_temp(isnan(uba_temp))=0;
ubo_temp(isnan(ubo_temp))=0;
meal_temp(isnan(meal_temp))=0;

uba = uba_temp;     % Filtered numerical values for basal insulin
ubo = ubo_temp;     % Filtered numerical values for bolus insulin
meal = meal_temp;   % Filtered numerical values for meals

%% Conversion factors

mmolLiter2mgDl = 18; % Multiply with 18 for converting mmol/L to mg/dL
h2min = 60;          % Convert from h to min
min2h = 1/h2min;     % Convert from min to h 
U2mU  = 1e3;         % Convert from U to mU 

%% Converting the measurements to correct units
CGM = CGM * mmolLiter2mgDl;    % Glucose measurements are given in mmol/L and we convert to mg/dL
uba = uba * U2mU*min2h;        % Insulin measurements are given in U/h and we convert to mU/min
ubo = ubo * U2mU*min2h;        % Insulin measurements are given in U/h and we convert to mU/min
%% Time format conversion

% We loop over the cell array and convert it to a string such that it gets
% the right format for datetime
for i = 1:length(time_temp_CGM)
    t(i) = datetime(time_temp_CGM{i},'InputFormat','dd/MM/yyyy HH:mm');
end

% Datetime for insulin and carb
for i = 1:length(time_temp_insulin_carb)
    t2(i)= datetime(time_temp_insulin_carb{i},'InputFormat','dd/MM/yyyy HH:mm');
end

%% Plots over the 4 measurements, G, ubo, uba and D for the whole time frame from january-june
% Plot glucose concentration from january till June
subplot(4,1,1)
plot(t, CGM)
ylim([30 400])
xlim([t(1) t(end)])
ylabel({'Glucose concentration', '[mg/gL]'});
xlabel('Time');
title('Blood glucose concentration over time')

% Plot basal insulin from january till june
subplot(4,1,2)
plot(t2, uba)
ylim([0 max(uba)])
xlim([t2(1) t2(end)])
ylabel({'Basal insulin', '[mU/min]'});
xlabel('Time');
title('Meals and meal sizes')

% Plot bolus insulin from january till june
subplot(4,1,3)
plot(t2, ubo)
ylim([0 max(ubo)])
xlim([t2(1) t2(end)])
ylabel({'Bolus insulin', '[mU/min]'});
xlabel('Time');
title('Basal insulin flow rate')

% Plot meals from january till june
subplot(4,1,4)
plot(t2, meal)
xlim([t2(1) t2(end)])
ylabel({'Meals', '[carbs in g]'});
xlabel('Time');
title('Bolus insulin')

%% Plots over the 4 measurements for almost a month
% We have approximately 130 days for t_datetime (length(t)/ 24*(60/5))
% approx 26 days a month, which is length(t)/5 = 7500 of 5-minute
% samplings. So we take the month march-april

t_test_CGM = t(floor((length(t)/5)*2):floor((length(t)/5)*2)+length(t)/5  ); % Time for almost month for CGM
t_test_ins_carb = t2(floor((length(t2)/5)*2):floor((length(t2)/5)*2)+length(t2)/5  ); % Time for almost month for insulin and carb

CGM_month = CGM(floor((length(CGM)/5)*2):floor((length(CGM)/5)*2)+length(CGM)/5);
uba_month = uba(floor((length(uba)/5)*2):floor((length(uba)/5)*2)+length(uba)/5);
ubo_month = ubo(floor((length(ubo)/5)*2):floor((length(ubo)/5)*2)+length(ubo)/5);
meal_month = meal(floor((length(meal)/5)*2):floor((length(meal)/5)*2)+length(meal)/5);

figure 
subplot(4,1,1)
plot(t_test_CGM, CGM_month)
ylim([min(CGM_month) max(CGM_month)])
xlim([t_test_CGM(1) t_test_CGM(end)])
ylabel({'CGM measurements', '[mg/dL]'});
xlabel('Time');
title('Blood glucose concentration over time')

% Plot meals from january till june
subplot(4,1,2)
plot(t_test_ins_carb, meal_month)
xlim([t_test_ins_carb(1) t_test_ins_carb(end)])
ylabel({'Meal carbohydrates', '[g CHO]'});
xlabel('Time');
title('Meals and meal sizes')

subplot(4,1,3)
plot(t_test_ins_carb, uba_month)
ylim([0 max(uba_month)])
xlim([t_test_ins_carb(1) t_test_ins_carb(end)])
ylabel({'Basal insulin', '[mU/min]'});
xlabel('Time');
title('Basal insulin flow rate')

% Plot bolus insulin from january till june
subplot(4,1,4)
plot(t_test_ins_carb, ubo_month)
ylim([0 max(ubo_month)])
xlim([t_test_ins_carb(1) t_test_ins_carb(end)])
ylabel({'Bolus insulin', '[mU]'});
xlabel('Time');
title('Bolus insulin')


%% GRID
delta_G        = 15;                 % From article
t_vec          = [5, 10,15];         % The respective sampling times
tau            = 6;                  % From the article
Gmin           = [90 0.5 0.5];       % For meal under 50 considered
Ts             = 5;                  % min - step size between control steps

D_detected = GRIDalgorithm_mealdetection(CGM_month,Gmin,tau,delta_G,t_vec,Ts);

% The total amount of detected meals
number_detectedmeals = sum(D_detected);

fprintf('number of detected meals: %d\n',number_detectedmeals);



