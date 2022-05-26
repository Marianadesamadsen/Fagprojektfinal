function Gf_ctrlstate = lowfilt_func(tau,tspan,Gfns,Gf_prev)
% 
% lowfilt_func()
% 
% DESCRIPTION:
% The function is a part of the GRID algortihm. This is part of the
% preprocessing section of the algortihm. The function filters the data 
% already filtered with the noise-spike filter. This means the function
% will damp the high frequency fluctuations. The function will smooth the
% data such that there will not be a long delay to optimize the specificity
% and detection speed of the algorithm.
% 
% INPUT:
% tau           -   Filter time constant
% delta_t       -   The samling time period 
% Gfns          -   The filtered value from spikefilt_func
% Gf_prev       -   The previous filtered value from lowfilt_func
%
% OUTPUT:
% The filtered value of the spikefilt_func value at the control state
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


% The low filter consists of one equation as follows
Gf_ctrlstate = ...
(tspan/(tau + tspan)) * Gfns + (1 - (tspan/ (tau + tspan)))* Gf_prev; 

end

