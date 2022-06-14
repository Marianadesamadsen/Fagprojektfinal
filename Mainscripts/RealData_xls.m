%% Loading real clinical data from excel where data is stored in sheets

clear all
clc

%% Adding the path 

addpath(genpath(fullfile(pwd, '../Other')));

CGM_sheet = readtable('Control-IQ_Sample_Diasend','Sheet','CGM');
Insulin_carbs_sheet = readtable('Control-IQ_Sample_Diasend', 'Sheet','Insulin and Carbs');

% Store the data in temporary arrays
CGM                     = table2array(CGM_sheet(:,2));              % Glucose concentration
uba_temp                = table2array(Insulin_carbs_sheet(:,2));    % Basal insulin
ubo_temp                = table2array(Insulin_carbs_sheet(:,4));    % Bolus insulin
meal_temp               = table2array(Insulin_carbs_sheet(:,8));    % Amount of carbs

% Store time in cell
time_temp_CGM           = table2cell(CGM_sheet(:,1));               % Time frame for glucose measurements
time_temp_insulin_carb  = table2cell(Insulin_carbs_sheet(:,1));     % Time frame for insulin and meal measurements

%% Replacing the missing values NaN with zeros

% Replacing all NaN values with 0
uba_temp(isnan(uba_temp))       = 0;
ubo_temp(isnan(ubo_temp))       = 0;
meal_temp(isnan(meal_temp))     = 0;

% Changning variables names after removing NAN with zeros
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

%% GRID - temp

% Preparing the time vector for GRID

TT = zeros(1,length(t));
TT(1) = 0;

% Making vector to numerical span between measurements
for i = 1:length(t)-1
    temp_time1 = datevec(t(i));
    temp_time2 = datevec(t(i+1));
    
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
Gmin           = [130 1.5 1.6];      % For meal under 50 considered

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
plot(t, CGM)
ylim([30 400])
xlim([t(1) t(end)])
ylabel({'Glucose concentration', '[mg/gL]'});
xlabel('Time');
hold on
plot(t(1:end-1),D_detected*100,'r.');

% Plot basal insulin from january till june
% subplot(4,1,2)
% plot(t2, uba)
% ylim([0 max(uba)])
% xlim([t2(1) t2(end)])
% ylabel({'Basal insulin', '[mU/min]'});
% xlabel('Time');
% 
% Plot bolus insulin from january till june
% subplot(4,1,3)
% plot(t2, ubo)
% ylim([0 max(ubo)])
% xlim([t2(1) t2(end)])
% ylabel({'Bolus insulin', '[mU/min]'});
% xlabel('Time');
% 
% Plot meals from january till june
% subplot(4,1,4)
% plot(t2, meal)
% xlim([t2(1) t2(end)])
% ylabel({'Meals', '[carbs in g]'});
% xlabel('Time');

% %% Plots over the 4 measurements for almost a month
% % We have approximately 130 days for t_datetime (length(t)/ 24*(60/5))
% % approx 26 days a month, which is length(t)/5 = 7500 of 5-minute
% % samplings. So we take the month march-april
% 
% t_test_CGM = t(floor((length(t)/5)*2):floor((length(t)/5)*2)+length(t)/5  ); % Time for almost month for CGM
% t_test_ins_carb = t2(floor((length(t2)/5)*2):floor((length(t2)/5)*2)+length(t2)/5  ); % Time for almost month for insulin and carb
% 
% CGM_month = CGM(floor((length(CGM)/5)*2):floor((length(CGM)/5)*2)+length(CGM)/5);
% uba_month = uba(floor((length(uba)/5)*2):floor((length(uba)/5)*2)+length(uba)/5);
% ubo_month = ubo(floor((length(ubo)/5)*2):floor((length(ubo)/5)*2)+length(ubo)/5);
% meal_month = meal(floor((length(meal)/5)*2):floor((length(meal)/5)*2)+length(meal)/5);
% 
% %% GRID
% delta_G        = 15;                 % From article
% tau            = 6;                  % From the article
% Gmin           = [130 1.5 1.6];      % For meal under 50 considered
% 
% % making vector to numerical span between measurements
% TT = zeros(1,length(t_test_CGM));
% TT(1) = 0;
% for i = 1:length(t_test_CGM)-1
%     temp_time1 = datevec(t_test_CGM(i));
%     temp_time2 = datevec(t_test_CGM(i+1));
%     
%     if temp_time2(4)==00 && temp_time1(4)~=0
%         TT(i+1) = TT(i)+(24*h2min+temp_time2(5)) - (temp_time1(4)*h2min+temp_time1(5));
%     else 
%         TT(i+1) = TT(i)+(temp_time2(4)*h2min+temp_time2(5)) - (temp_time1(4)*h2min+temp_time1(5));
%     end
% end
% 
% % step size between control steps
% Ts = zeros(1,length(TT));
% for k = 1:length(TT)-1
%     Ts(k) = TT(k+1)-TT(k);
% end
% 
% D_detected = GRIDalgorithm_mealdetection2(CGM_month,Gmin,tau,delta_G,TT,Ts);
% 
% % The total amount of detected meals 
% number_detectedmeals = sum(D_detected);
% 
% fprintf('number of detected meals: %d\n',number_detectedmeals);
% 
% %% Visualize 
% figure;
% 
% subplot(4,1,1)
% plot(t_test_CGM, CGM_month)
% ylim([min(CGM_month) max(CGM_month)])
% xlim([t_test_CGM(1) t_test_CGM(end)])
% ylabel({'Glucose concentration', '[mg/dL]'});
% xlabel('Time');
% hold on
% plot(t_test_CGM(1:end-1),D_detected*100,'r.');
% 
% subplot(4,1,2)
% stairs(t_test_ins_carb, uba_month)
% ylim([0 max(uba_month)])
% xlim([t_test_ins_carb(1) t_test_ins_carb(end)])
% ylabel({'Basal insulin', '[mU/min]'});
% xlabel('Time');
% 
% % Plot bolus insulin from january till june
% subplot(4,1,3)
% stem(t_test_ins_carb, ubo_month*mU2U, 'markersize',1)
% %ylim([0 max(ubo_month)])
% xlim([t_test_ins_carb(1) t_test_ins_carb(end)])
% ylabel({'Bolus insulin', '[U/min]'});
% xlabel('Time');
% 
% % Plot meals from january till june
% subplot(4,1,4)
% plot(t_test_ins_carb, meal_month*Ts)
% xlim([t_test_ins_carb(1) t_test_ins_carb(end)])
% ylabel({'Meals', '[g CHO]'});
% xlabel('Time');
% 
% 
% 
