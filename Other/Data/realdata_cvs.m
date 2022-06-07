%% Computing GRID algorithm on real data 

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

%% Loading data

data=importdata('Control-IQ_Sample_Tconnect.csv');

% Extracting the glucose concentration
G = data.data;
% Extracting the time
date_temp = data.textdata(2:780,4);

% Removing the T's in the dates and replaces with a space such that it gets the right format in datetime
date = regexprep(date_temp, 'T', ' '); 

% We loop over the cell array and convert it to a string such that it gets
% the right format for datetime
for i = 1:length(date)
    str = string(date{i});
    t(i)= datetime(str,'InputFormat','yyyy-MM-dd HH:mm:ss');
end

%% Visualize 

% Plot over real patient's glucose levels over 3 days
figure (1)
subplot(2,1,1)
plot(t, G)
title('Clinical patient glucose conc. over 3 days')
xlabel('Time')
ylabel('Blood glucose concentration')

%% Detecting meals using GRID algorithm

% Inisializing
delta_G        = 15;                 % From article
t_vec          = [5, 10,15];         % The respective sampling times
tau            = 6;                  % From the article
Gmin           = [90 0.5 0.5];       % For meal under 50 considered
Ts             = 5;                  % min - step size between control steps

D_detected = GRIDalgorithm_mealdetection(G,Gmin,tau,delta_G,t_vec,Ts);

% The total amount of detected meals
number_detectedmeals = sum(D_detected);

fprintf('number of detected meals: %d\n',number_detectedmeals);

%% Visualize

% Plot blood glucose concentration and the detected meals as points
subplot(2,1,2)
plot(t, G);
ylabel({'Blood glucose concentration', '[mg/dL]'});
xlabel('time');
hold on
plot(t(1:end-1),D_detected*150,'r.');
title('Clinical patient glucose conc. over 3 days with detected meals in red dots')


