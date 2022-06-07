function Gfm_ctrlstate = estimate_lagrange(t_vec,Gf_vec)
% 
% estimate_lagrange()
% 
% DESCRIPTION:
% The function is a part of the GRID algortihm. This is part of the
% estimation section of the algortihm. The function finds the first
% derivative using the 3-point lagrangian interpolation polynomial. The
% first derivative is used later in the last section of the GRID algorithm.
% 
% INPUT:
% t_vec   - the sampling time at t(k),t(k-1),t(k-2) in a vector
% Gf_vec  - the lowfilt_func data at the above time samples in a vector
%
% OUTPUT:
% The first derivative of the filtered glucose measurement at the control
% sate
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

% Inizialising the values of the sampling times 
tk      = t_vec(3);
tkm1    = t_vec(2);
tkm2    = t_vec(1);

% Inizialising the values of the filtered glucose measmurements
Gf      = Gf_vec(3);
Gfm1    = Gf_vec(2);
Gfm2    = Gf_vec(1);

% Calculating the first derivative at t(k)
Gfm_ctrlstate = ...
      ( (tk - tkm1)/( (tkm2-tkm1)*(tkm2-tk) ) ) * Gfm2 ...
    + ( (tk - tkm2)/( (tkm1-tkm2)*(tkm1-tk) ) ) * Gfm1 ...
    + ( (2*tk - tkm2 - tkm1)/( (tk-tkm1)*(tk-tkm2) ) ) * Gf;

end
