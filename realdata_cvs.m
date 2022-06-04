%% Simulation test of the GRID algorithm 30 days, 3 meals with noise

% Simulating 3 meals on 30 days, and detecting the meals
% by using the GRID algorithm

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
clc 
clear

data=importdata('Control-IQ_Sample_Tconnect.csv');

G = data.data;
date_temp = data.textdata(2:780,4);

date = regexprep(date_temp, 'T', ' '); % Removes the T's in the dates and replaces with a space such that it gets the right format in datetime

% We loop over the cell array and convert it to a string such that it gets
% the right format for datetime
for i = 1:length(date)
    str = string(date{i});  
    t(i)= datetime(str,'InputFormat','yyyy-MM-dd HH:mm:ss'); 
end

% Plot over real patient's glucose levels over 3 days
figure
plot(t, G)
title('Clinical patient glucose conc. over 3 days')
xlabel('Time')
ylabel('Blood glucose concentration')

%% Number of control steps
Ts = 5;         % min - step size 

%% Detecting meals using GRID algorithm

% Inisializing 
filt_prev      = zeros(length(G),2); % The vector of previous filtered values
Gfm_vec        = zeros(length(G),2); % The vector of previous derivatives
G_vec          = [G(1),G(1),G(1)];   % Inserting the start previous glucose measurements as the same value
delta_G        = 15;                 % From article
t_vec          = [5, 10,15];   % The respective sampling times
filt_prev(1,:) = [G(1),G(1)];        % Inserting the previous filtered value as the not filtered values
tau            = 6;                  % From the article
flag           = 0;                  % No detected meal to begin with
Gmin           = [120 0.8 0.65];     % Gmin accepts intensity up to 6

% Other tries
%Gmin           = [90 0.5 0.5];       % For meal under 50 considered
%Gmin = [ 130 1.5 1.6 ]; % Their meals
%Gmin = [ 110 1 1.5 ]; % For no meal under 50 


% Computing the first two detections
[ Gfm_vec(2,:) , filt_prev(2,:) , flag, zero_one(2) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(1,:) , Gmin, Gfm_vec(1,:) , t_vec ,flag );
                                
                                % Updating G_vec
                                G_vec=[G(1),G(1),G(2)];
                                
[ Gfm_vec(3,:) , filt_prev(3,:) , flag, zero_one(3) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(2,:) , Gmin, Gfm_vec(2,:) , t_vec, flag );
                                
                                % Updating G_vec
                                G_vec=[G(1),G(2),G(3)];
                                
% Computing the last detections
for i = 3 : length(G)-1
    
[ Gfm_vec(i+1,:) , filt_prev(i+1,:) , flag, zero_one(i) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(i,:) , Gmin, Gfm_vec(i,:) , t_vec , flag );
                                
                                % Updating G_vec
                                G_vec=[G(i-1),G(i),G(i+1)];
                                
end

% The total amount of detected meals
detectedmeals=sum(zero_one);

%% Visualize 

% Create figure with absolute size for reproducibility
figure;

% Plot blood glucose concentration and the detected meals as points
subplot(211);
plot(t, G);
ylabel({'Blood glucose concentration', '[mg/dL]'});
hold on 
plot(t(1:end-1),zero_one*150,'r.');

% Plot detected Meals
subplot(212);
plot(t(1:end-1),zero_one,'b-');
ylabel({'detected meal'}); 
xlabel('Time [h]'); 
 

