%%

csvread('Control-IQ_Sample_Tconnect.csv', 2,2)

opts.SelectedVariableNames = [1:5];

opts.Sheet = 'CGM'
readmatrix('Control-IQ_Sample_Tconnect.csv')

%%
data_temp = readtable('Control-IQ_Sample_Tconnect.csv');

data = data_temp(1:779, 3);

%% 
data=importdata('Control-IQ_Sample_Tconnect.csv')


