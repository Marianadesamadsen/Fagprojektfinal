%%

csvread('Control-IQ_Sample_Tconnect.csv', 2,2)

opts.SelectedVariableNames = [1:5];

opts.Sheet = 'CGM'
readmatrix('Control-IQ_Sample_Tconnect.csv')

%%
data_temp = readtable('Control-IQ_Sample_Tconnect.csv');

data = data_temp(1:779, 3);

%% 
% Load the full data
data_temp = importdata('Control-IQ_Sample_Tconnect.csv');
data = struct2cell(data_temp); 

% Convert to CGM data and date of time
CGM_data = data{1}

date_temp = data{2}
date = date_temp(1:780, 4)  








