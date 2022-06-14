function pmat = pmatrix(n)
%
% pmatrix()
% 
% DESCRIPTION:
% Simulates a P matrix over n number of patients. 
% 
% INPUT: 
% p         - Number of patients 
% 
% OUTPUT:
% Matrix with p vectors for every number of patients. Dim: p x n.
%
% PROJECT:
% Fagprojekt 2022
% A diabetes case study - Meal detection
%
% GENEREL:
% BSc                       : Mathematics and technology 
% University                : The Technical University of Denmark (DTU)
% Department                : Applied Mathematics and Computer Science 
% 
% AUTHORS:
% Emma Victoria Lind
% Mariana de Sá Madsen 
% Mona Saleem
% 
% CONTACT INFORMATION
% s201205@student.dtu.dk
% s191159@student.dtu.dk
% s204226@student.dtu.dk
% 
% REFERENCE:
% MANGLER FRA ARTIKEL
%

% The different parameters being simulated in range from article 
p1 = 41+(131-41).*rand(1,n);
p2 = 10+(70-10).*rand(1,n);
p3 = 5.40+(20.10-5.40).*rand(1,n);
p4 = 0.0081+(0.0233-0.0081).*rand(1,n);
p5 = (9.64e-4+(1.7e-2-9.64e-4).*rand(1,n));
p6 = 1*10^(-8)+(6.39*10^(-3)-1*10^(-8)).*rand(1,n);
p7 = 0.6+(3.45-0.6).*rand(1,n);
p8 = 104+(337-104).*rand(1,n);
p9 = 47*ones(1,n);
p10 = 5*ones(1,n);

% Every vector set up in matrix
pmat = [p1;p2;p3;p4;p5;p6;p7;p8;p9;p10];

end






