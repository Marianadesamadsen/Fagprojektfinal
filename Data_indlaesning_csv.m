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
CGM_data = data{1};

date_temp = data{2};
date = date_temp(2:780, 4);  % The date in cell aray


DATE = regexprep(date, 'T', ' '); % Removes the T's in the dates and replaces with a space such that it gets the right format in datetime

% We loop over the cell array and convert it to a string such that it gets
% the right format for datetime
for i = 1:length(DATE)
    str = string(DATE{i});  
    t(i)= datetime(str,'InputFormat','yyyy-MM-dd HH:mm:ss');
end

% Test med forskellige 'funktioner'
% date_string = cellstr(date); % The date in cell aray
%date_char = cell2mat(date) % The date in character array

% for i = 1:length(DATE)
%     str = 
% d = '2018-06-25 11:23:37.712';
% t = datetime(test,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS')


% %% EVENTUELT BRUG DEN HER NEDENSTÅENDE FOR DER BEHØVES MAN IKKE FJERNE T
% !!!!!!
% t = datetime(date_char,'InputFormat','uuuu-MM-dd HH:mm:ss''T''HH:mmXXX','TimeZone','UTC')

%%
XDates = [datetime(2021,6,1:30) datetime(2021,7,1:31)];
YNumsForXDates = sin(rand(1,length(t)));
plot(t,YNumsForXDates)
