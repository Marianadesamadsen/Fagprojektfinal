%% Loading real clinical data from excel where data is stored in sheets
clear all
clc
close all
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
CGM = CGM*mmolLiter2mgDl;    % Glucose measurements are given in mmol/L and we convert to mg/dL
uba = uba*U2mU*min2h;        % Insulin measurements are given in U/h and we convert to mU/min
ubo = ubo*U2mU*min2h;        % Insulin measurements are given in U/h and we convert to mU/min
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

%% Plots over approx 3 days
% Plot glucose concentration from january till June
plot(t(1:720), CGM(1:720))
ylim([30 400])
xlim([t(1) t(720)])
ylabel({'Glucose concentration', '[mg/dL]'});
xlabel('Time');
