function g = CGMsensor(x,pg)
%
% CGMsensor()
%
% DESCRIPTION:
% Subtracts the glucose value from the x vector
%
% INPUT:
% x     - the state vector
% pg    - vector of parameters (not used)
%
% OUTPUT:
% The glucose value
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

% [mg/dL] The glucose concentration
g = x(6,:);

end
