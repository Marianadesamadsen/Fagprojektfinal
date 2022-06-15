%% Plot over 3 days from january 2nd till end of January 4th
t_start2 = t2(9734)    % Start time, the one out of the two measurements that start first
t_end2 = t2(10227) % Finish time, the one out of the two measurements that finishes last

time_vec_CGM2 = [t_start2, t(16297:17156),t_end2]         % Datetime vector for CGM (add one more timepoint)
time_vec_insulin2 = t2                     % Datetime vector for insulin & carbs (add one more timepoint)

time3days_CGM = time_vec_CGM2;                     % Datetime vector for 3 days for CGM
time3days_insulin = time_vec_insulin2(9734:10227) % Datetime vector for 3 days for insulin & carbs

CGM_32 = CGM(16297:17156)       % CGM measurements for 3 days
CGM_32(end+1:end+2) = 0         % Add an extra row with 0 to match the dimension of timevector

ubo_32 = ubo(9734:10227)       % Bolus measurements for 3 days

uba_32 = uba(9734:10227)       % Basal measurements for 3 days

meal_32 = meal(9734:10227)     % Meal measurements for 3 days

%% GRID over 3 days
TT = zeros(1,length(time3days_CGM));
TT(1) = 0;

% Making vector to numerical span between measurements
for i = 1:length(time3days_CGM)-1
    temp_time1 = datevec(time3days_CGM(i));
    temp_time2 = datevec(time3days_CGM(i+1));

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
Gmin           = [120 0.6 0.5];      % For meal under 50 considered

% Computing the GRID algorithm
D_detected = GRIDalgorithm_mealdetection2(CGM_32,Gmin,tau,delta_G,TT,Ts);

% The total amount of detected meals
number_detectedmeals = sum(D_detected);

% Printing the number of detected meals
fprintf('number of detected meals: %d\n',number_detectedmeals);

D_detected(end+1)=0; % Adding one extra row to match the dimension of time vector


%% Plots over 3 days
figure (2)
subplot(4,1,1)
plot(time3days_CGM, CGM_32);
hold on
plot(time3days_CGM, D_detected*200, '*')
xlim([t(16297) t(17156)])
ylabel({'Blood glucose concentration', '[mg/dL]'});
xlabel('Time');
title('Blood glucose concentration over time')

subplot(4,1,2)
stem(time3days_insulin, ubo_32*h2min*mU2U)
ylabel({'Bolus insulin', '[U]'});
xlabel('Time');
title('Bolus')

subplot(4,1,3)
stairs(time3days_insulin, uba_32)
ylabel({'Basal insulin', '[mU/min]'});
xlabel('Time');
title('Basal insulin flow rate')

subplot(4,1,4)
stem(time3days_insulin, meal_32)
ylabel({'Meal carbohydrates', '[g CHO]'});
xlabel('Time');
title('Meals and meal sizes')