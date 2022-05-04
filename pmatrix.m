function pmat = pmatrix(numberpatients)
% P parameter matrix function 
% 
% GOAL:
% Simulates a P matrix over number patients. dim  = p x numberpatients
% 
% col: patients
% row: parameter i
% 
% INPUT: 
% Number of patients desired to be simulated
% 
% OUTPUT:
% Matrix with p vectors for every number of patients


% The different parameters being simulated in range from article 
p1 = 41+(131-41).*rand(1,numberpatients);
p2 = 10+(70-10).*rand(1,numberpatients);
p3 = 540+(2010-540).*rand(1,numberpatients);
p4 = 0.0081+(0.0233-0.0081).*rand(1,numberpatients);
p5 = 9.6400e-05+(9.48*10^(-4)-9.6400e-05).*rand(1,numberpatients);
p6 = 1*10^(-8)+(6.39*10^(-3)-1*10^(-8)).*rand(1,numberpatients);
p7 = 0.6+(3.45-0.6).*rand(1,numberpatients);
p8 = 104+(337-104).*rand(1,numberpatients);
p9 = 47*ones(1,numberpatients);
p10 = 5*ones(1,numberpatients);

% Every vector set up in matrix
pmat = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10];

end






