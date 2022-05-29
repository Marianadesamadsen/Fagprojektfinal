function h = CHMsensor_withnoise(X,ph)
% 
% CGMsensor()
% 
% DESCRIPTION:
% Subtracts the subcutaneous glucose value from the x vector and adds gaussian noise
%
% INPUT:
% x     - the state vector
% ph    - vector of parameters (not used)
%
% OUTPUT:
% h     - the subcutaneous blood glucose concentration with Gaussian noise
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

% Subcutanius glucose concentration
h = X(7, :);

% Computing the Gaussian noise
std=1.5;        % standard deviation of 2%
meanValue=0;    % mean=0

% Adding the noise
Noise_signal_g = h + std*randn(size(h)) + meanValue;

% Updating output with noise
h = Noise_signal_g;

end



