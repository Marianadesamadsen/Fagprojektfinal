%%

csvread('Control-IQ_Sample_Tconnect.csv', 2,2)

opts.SelectedVariableNames = [1:5];
readmatrix('Control-IQ_Sample_Tconnect.csv')

%%
data = readtable('Control-IQ_Sample_Tconnect.csv')
