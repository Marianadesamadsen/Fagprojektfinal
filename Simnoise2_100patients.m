%% Simulation test of the GRID algorithm on 100 patients with stochastic measurement noise

% Simulating 3 meals and 2 snacks on 30 days on 100 patients at a time
% using the GRID algorithm with updated stochastic measurement noise

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
min2sec = h2min; % Convert from min to sec
sec2min = 1/h2min;% Convert from sec to min

%% Inizializing parameters

numpatients = 100;         % number of patients
pf = pmatrix(numpatients); % computing the p vectors for all patients
Gs = 108;                  % Steady state blood glucose concentration
ts = [];

%% Computing steadty state

% Initializing
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

[T(:,:,p), X(:,:,p)] = OpenLoopSimulation(x0(p,:)', tspan, U(:,:,p), D(:,:,p), pf(:,p), @MVPmodel, @EulerM, Nk);

end 


%% Blood glucose concentration 

% Inisializing
G=zeros(1,N+1,numpatients);

% Looping over all patients
for p=1:numpatients

G(:,:,p) = CGMsensor_withnoise(X(:,:,p), pf(:,p)); % [mg/dL] 

end

%% GRID

% Inisializing
filt_prev      = zeros(length(G(:,:,1)),2,numpatients); % The vector of previous filtered values
Gfm_vec        = zeros(length(G),2,numpatients);        % The vector of previous derivatives
detectedmeals  = zeros(1,numpatients);                  % Detected meals for each patient in matrix
zero_one       = zeros(length(G)-1,numpatients);        % Detected meals and not detected meals 

for p = 1:numpatients 

G_vec            = [G(:,1,p),G(:,1,p),G(:,1,p)]; % Inserting the start previous glucose measurements as the same value
delta_G          = 15;                           % From article  
t_vec            = [5,10,15];                    % The respective sampling times
filt_prev(1,:,p) = [G(:,1,p),G(:,1,p)];          % Inserting the previous filtered value as the not filtered values
tau              = 6;                            % From article
flag             = 0;                            % No detected meals to begin with 
Gmin             = [90 0.5 0.5];                 % For meal under 50 considered

% Other tries
% Gmin = [ 130 1.5 1.6 ]; % Their meals 
% Gmin = [ 110 1 1.5 ]; % For no meal under 50 

% Computing the first two detections
[ Gfm_vec(2,:,p) , filt_prev(2,:,p) , flag, zero_one(1,p) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(1,:,p) , Gmin, Gfm_vec(1,:,p) , t_vec ,flag );
                                
                                % Updating G_vec
                                G_vec=[G(:,1,p),G(:,1,p),G(:,2,p)];
                                
[ Gfm_vec(3,:,p) , filt_prev(3,:,p) , flag, zero_one(2,p) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(2,:,p) , Gmin, Gfm_vec(2,:,p) , t_vec, flag );
                                
                                % Updating G_vec
                                G_vec=[G(:,1,p),G(:,2,p),G(:,3,p)];

        % Comptuing the last detections                               
        for i = 3 : length(G)-1

         [ Gfm_vec(i+1,:,p) , filt_prev(i+1,:,p) , flag, zero_one(i,p) ]= ...
             GRID_func( delta_G , G_vec , tau, Ts , ...
             filt_prev(i,:,p) , Gmin, Gfm_vec(i,:,p) , t_vec , flag );
                   
         % Updating G_vec
         G_vec=[G(i-1),G(i),G(i+1)];
                                
        end
         
% The total amount of detected meals for each patient in vector
detectedmeals(p)=sum(zero_one(:,p));

end

%% Finding minimum and maximum patient

% Sums the glucose concentration for each patient and stores it as a vector
s=sum(G(:,:,1:100));

% Find the minimum value
minp=min(s);

% Find the maximum value
maxp=max(s);

% Loops over alle possible sums to find the index for the maximum and
% minimum patient
for i=1:100
    
    % If the ith patient has the minimum patient sum then it is the minum
    % patient
    if minp == sum(G(:,:,i))
        minpatient = i;
        minpatient = i; % The index for the max patient
    end
    
    % If the ith patient has the maximum patient sum then it is the maximum
    % patient
    if maxp == sum(G(:,:,i))
        maxpatient = i; % The index for the max patient
    end
    
end

%% Visualize 

reset(groot);

% Create figure with absolute size for reproducibility
figure;

% Converting data
T2=datetime(T*min2sec,'ConvertFrom','posixtime');
tspan2=datetime(tspan*min2sec,'ConvertFrom','posixtime');

% Plot blood glucose concentration and the detected meals as points
subplot(411);
plot((T2(:,:,minpatient)), G(:,:,minpatient),'r');
hold on 
plot((T2(:,:,maxpatient)), G(:,:,maxpatient),'b');
%xlim([t0, tf]*min2h);
ylabel({'Blood glucose concentration', '[mg/dL]'});
hold on 
plot(tspan2(1:end-1),zero_one(:,minpatient)*400,'*r');
hold on 
plot(tspan2(1:end-1),zero_one(:,maxpatient)*300,'*b');
legend('minpatient','maxpatient','minpatient','maxpatient');

% Plot meal carbohydrate and the detected meals as points
subplot(412);
stemmin = Ts*D(1,:,minpatient);
stemmax = Ts*D(1,:,maxpatient);
stem(tspan2(1:end-1),stemmin,stemmax, 'MarkerSize', 0.1);
% hold on 
% stem(tspan2(1:end-1), Ts*D(1,:,maxpatient), 'MarkerSize', 0.1);
%xlim([t0, tf]*min2h);
ylabel({'Meal carbohydrates', '[g CHO]'});
hold on 
plot(tspan2(1:end-1),zero_one(:,minpatient)*150,'*r');
hold on 
plot(tspan2(1:end-1),zero_one(:,maxpatient)*100,'*b');
legend('minpatient','maxpatient','minpatient','maxpatient');

% Plot basal insulin flow rate
subplot(413);
stairs(tspan2, U(1, [1:end, end],minpatient));
stairs(tspan2, U(1, [1:end, end],maxpatient));
legend('minpatient','maxpatient');
%xlim([t0, tf]*min2h);
ylabel({'Basal insulin', '[mU/min]'});

% Plot bolus insulin
subplot(414);
stem(tspan2(1:end-1), Ts*mU2U*U(2, :,minpatient), 'MarkerSize', 1);
stem(tspan2(1:end-1), Ts*mU2U*U(2, :,maxpatient), 'MarkerSize', 1);
%xlim([t0, tf]*min2h);
ylabel({'Bolus insulin', '[U]'}); 
xlabel('Time [h]');










