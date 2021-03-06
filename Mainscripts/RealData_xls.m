%% Real clinical data
% Loading of real clinical data over 5 months given by an excel file 
% where data is  stored in sheets.
% Extracting the relevant measurements: Glucose concenctration by CGM,
% bolus insulin, basal insulin and size of meals given in gram carbs and
% the time for each measurement. The measurements are plotted and the
% GRID-algorithm is tested on the real data.

%%
clear all
clc

%% Adding the path and loading the file and measurements

addpath(genpath(fullfile(pwd, '../Other')));

CGM_sheet = readtable('Control-IQ_Sample_Diasend','Sheet','CGM');
Insulin_carbs_sheet = readtable('Control-IQ_Sample_Diasend', 'Sheet','Insulin and Carbs');

% Store the data in temporary arrays
CGM                     = table2array(CGM_sheet(:,2));              % Glucose concentration
uba_temp                = table2array(Insulin_carbs_sheet(:,2));    % Basal insulin
ubo_temp                = table2array(Insulin_carbs_sheet(:,4));    % Bolus insulin
meal_temp               = table2array(Insulin_carbs_sheet(:,8));    % Amount of carbs

% Store time in temporary cells
time_temp_CGM           = table2cell(CGM_sheet(:,1));               % Time frame for glucose measurements
time_temp_insulin_carb  = table2cell(Insulin_carbs_sheet(:,1));     % Time frame for insulin and meal measurements

%% Filtering the data by replacing the missing values NaN with zeros

% Replacing all NaN values with 0
uba_temp(isnan(uba_temp))       = 0;
ubo_temp(isnan(ubo_temp))       = 0;
meal_temp(isnan(meal_temp))     = 0;

% Changning variables names after removing NaN with zeros
uba     = uba_temp;             % Filtered numerical values for basal insulin
ubo     = ubo_temp;             % Filtered numerical values for bolus insulin
meal    = meal_temp;            % Filtered numerical values for meals

%% Conversion factors

mmolLiter2mgDl  = 18;          % Multiply with 18 for converting mmol/L to mg/dL
h2min           = 60;          % Convert from h to min
min2h           = 1/h2min;     % Convert from min to h
U2mU            = 1e3;         % Convert from U to mU
mU2U            = 1/U2mU;      % Convert from U to mU

%% Converting the measurements to correct units

CGM = CGM * mmolLiter2mgDl;    % Glucose measurements are given in mmol/L and we convert to mg/dL
uba = uba * U2mU*min2h;        % Insulin measurements are given in U/h and we convert to mU/min
ubo = ubo * U2mU*min2h;        % Insulin measurements are given in U/h and we convert to mU/min

%% Time format conversion to datetime

% We loop over the cell array and convert it to a string such that it gets
% the right format for datetime

% Datetime for CGM measurements
for i = 1:length(time_temp_CGM)
    t(i) = datetime(time_temp_CGM{i},'InputFormat','dd/MM/yyyy HH:mm');
end

% Datetime for insulin and carb measurements
for i = 1:length(time_temp_insulin_carb)
    t2(i)= datetime(time_temp_insulin_carb{i},'InputFormat','dd/MM/yyyy HH:mm');

end

%% Data from january 2nd till end of January 4th
t_start = min(t(1),t2(1));    % Start time, the one out of the two measurements that start first
t_end = max(t(end), t2(end)); % Finish time, the one out of the two measurements that finishes last

time_vec_CGM = [t_start, t, t_end];         % Datetime vector for CGM (add one more timepoint)
time_vec_insulin = t2;                      % Datetime vector for insulin

time3days_CGM = time_vec_CGM(1:858+1);      % Datetime vector for 3 days for CGM
time3days_insulin = time_vec_insulin(1:482); % Datetime vector for 3 days for insulin & carbs

CGM_3 = CGM(1:858);       % CGM measurements for 3 days
CGM_3(end+1) = 0;         % Add an extra row with 0 to match the dimension of timevector

ubo_3 = ubo(1:482);       % Bolus measurements for 3 days

uba_3 = uba(1:482);       % Basal measurements for 3 days

meal_3 = meal(1:482);     % Meal measurements for 3 days

%% GRID on clinical data of 3 days
TT = zeros(1,length(time3days_CGM));
TT(1) = 0;

% Making vector to numerical span between measurements
for i = 1:length(time3days_CGM)-1
    temp_time1 = datevec(time3days_CGM(i));   % Convert datetime to vector of components with dates and hours
    temp_time2 = datevec(time3days_CGM(i+1)); % Convert datetime to vector of components with dates and hours

    if temp_time2(4)==00 && temp_time1(4)~=0
        TT(i+1) = TT(i)+(24*h2min+temp_time2(5)) - (temp_time1(4)*h2min+temp_time1(5));
    else
        TT(i+1) = TT(i)+(temp_time2(4)*h2min+temp_time2(5)) - (temp_time1(4)*h2min+temp_time1(5));
    end
end

% step size between control steps
Ts = zeros(1,length(TT));

for k = 1:length(TT)-1
    Ts(k) = TT(k+1)-TT(k); % This is the final vector to use
end

% Initializing for the GRID algorithm
delta_G        = 15;                 % From article
tau            = 6;                  % From the article
Gmin           = [120 0.6 0.5];      % The optimal found Gmin values

% Computing the GRID algorithm
D_detected = GRIDalgorithm_mealdetection2(CGM_3,Gmin,tau,delta_G,TT,Ts);

% The total amount of detected meals
number_detectedmeals = sum(D_detected);

% Printing the number of detected meals
fprintf('number of detected meals: %d\n',number_detectedmeals);

D_detected(end+1)=0; % Adding one extra row to match the dimension of time vector

figure(3)

subplot(4,1,1)
plot(time3days_CGM, CGM_3);
hold on
plot(time3days_CGM, D_detected*200, '*')
ylabel({'Blood glucose concentration', '[mg/dL]'});
title('Blood glucose concentration over time', 'FontSize', 25)

subplot(4,1,4)
stem(time3days_insulin, meal_3)
ylabel({'Meal carbohydrates', '[g CHO]'});
xlabel('Time [h]');
title('Meals and meal sizes', 'FontSize', 25)

subplot(4,1,2)
stem(time3days_insulin, ubo_3*h2min*mU2U)
ylabel({'Bolus insulin', '[U]'});
title('Bolus insulin', 'FontSize', 25)

subplot(4,1,3)
stairs(time3days_insulin, uba_3)
ylabel({'Basal insulin', '[mU/min]'});
title('Basal insulin flow rate', 'FontSize', 25)

%% Most optimal Gmin values tried on the GRID algorithm (several combinations)

% Allocate space for detections for each G_min
detec = zeros(length(idx_optimalfinal),length(CGM_3)-1);
detecmeals = zeros(length(idx_optimalfinal), 1);

% Define the table of optimal G min combinations and the indices
delta_G        = 15;                 % From article
tau            = 6;                  % From the article

% Total sum of detected meals for the most optimal Gmin values found from
% the simulated 100 patients
for i = 1:length(idx_optimalfinal)
    G_min_temp = Gmin_combinations(idx_optimalfinal(i),:);
    detec(i,:) = GRIDalgorithm_mealdetection2(CGM_3,G_min_temp,tau,delta_G,TT,Ts);
    detecmeals(i) = sum(detec(i,:));
end

detectedmeals_vec = detec(1,:);
detectedmeals = sum(detectedmeals_vec);

%% Plots over the 4 measurements from 1st of March till 31st of March

figure(3)

%Plot glucose concentration from 1st of march till 31st of march
subplot(4,1,1)
plot(t(11769:20536), CGM(11769:20536))
xlim([t(11767) t(20555)])
ylabel({'Glucose concentration', '[mg/gL]'});
title('Blood glucose concentration over time', 'FontSize', 25)


%Plot meals from 1st of march till 31st of march
subplot(4,1,2)
plot(t2(6990:12407), meal(6990:12407))
xlim([t2(6980) t2(12415)])
ylabel({'Meal carbohydrates', '[g CHO]'});
title('Meals and meal sizes', 'FontSize', 25)

%Plot basal insulin from 1st of march till 31st of march
subplot(4,1,3)
plot(t2(6990:12407), uba(6990:12407))
ylim([0 max(uba)])
xlim([t2(6980) t2(12415)])
ylabel({'Basal insulin', '[mU/min]'});
title('Basal insulin flow rate', 'FontSize', 25)

%Plot bolus insulin from 1st of march till 31st of march
subplot(4,1,4)
plot(t2(6990:12407), ubo(6990:12407)*h2min*mU2U)
xlim([t2(6980) t2(12415)])
ylabel({'Bolus insulin', '[U]'});
xlabel({'Time [h]'})
title('Bolus insulin', 'FontSize', 25)
    
%% ***** ALL OF THE FOLLOWING IS FOR THE WHOLE TIME FRAME FROM JANUARY 1ST TILL JUNE *****

%% GRID on clinical data for the whole timeframe of 5 months

% Preparing the time vector for GRID
TT = zeros(1,length(t));
TT(1) = 0;

% Making vector to numerical span between measurements
for i = 1:length(t)-1
    temp_time1 = datevec(t(i));   % Convert datetime to vector of components with dates and hours
    temp_time2 = datevec(t(i+1)); % Convert datetime to vector of components with dates and hours

    if temp_time2(4)==00 && temp_time1(4)~=0
        TT(i+1) = TT(i)+(24*h2min+temp_time2(5)) - (temp_time1(4)*h2min+temp_time1(5));
    else
        TT(i+1) = TT(i)+(temp_time2(4)*h2min+temp_time2(5)) - (temp_time1(4)*h2min+temp_time1(5));
    end
end

% step size between control steps
Ts = zeros(1,length(TT));

for k = 1:length(TT)-1
    Ts(k) = TT(k+1)-TT(k); % This is the final vector to use
end

% Initializing for the GRID algorithm
delta_G        = 15;                 % From article
tau            = 6;                  % From the article
Gmin           = [120 0.6 0.5];      % Optimal found Gmin values

% Computing the GRID algorithm
D_detected = GRIDalgorithm_mealdetection2(CGM,Gmin,tau,delta_G,TT,Ts);

% The total amount of detected meals
number_detectedmeals = sum(D_detected);

% Printing the number of detected meals
fprintf('number of detected meals: %d\n',number_detectedmeals);

%% Plots over the 4 measurements, G, ubo, uba and D for the whole time frame from january-june

%Plot glucose concentration from january till June
%subplot(4,1,1)
figure;
plot(t(1:7200), CGM(1:7200))
ylim([30 400])
xlim([t(1) t(end)])
ylabel({'Glucose concentration', '[mg/gL]'});
xlabel('Time');

Plot basal insulin from january till june
subplot(4,1,2)
plot(t2, uba)
ylim([0 max(uba)])
xlim([t2(1) t2(end)])
ylabel({'Basal insulin', '[mU/min]'});
xlabel('Time');

Plot bolus insulin from january till june
subplot(4,1,3)
plot(t2, ubo*mU2U)
ylim([0 max(ubo)])
xlim([t2(1) t2(end)])
ylabel({'Bolus insulin', '[U/min]'});
xlabel('Time');

Plot meals from january till june
subplot(4,1,4)
plot(t2, meal)
xlim([t2(1) t2(end)])
ylabel({'Meals', '[carbs in g]'});
xlabel('Time');


