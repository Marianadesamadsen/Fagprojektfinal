function h = CHMsensor(x,ph)
% CHMsensor()
% 
% DESCRIPTION:
% Function CHMsensor returns the value for the blood glucose concentration
% G(t) from the statevector.
%
% INPUT:
% x     - state vector x  (dimension: 7)
% ph    - parametervalues (dimension: 10)
%
% OUTPUT:
% Blood glucose concentration G(t)
% 
% PROJECT:
% Fagprojekt 2022
% A diabetes case study - Meal detection
%
% GENERAL:
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

% [mg/dL] Blood glucose concentration
    h = x(6);
    
end



