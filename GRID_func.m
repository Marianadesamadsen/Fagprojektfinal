function  [ Gfm_vec , G_prev , flag, zero_one ] = GRID_func( ...
         delta_G , G , tau, tspan , G_prev , Gmin, Gfm_vec , t_vec, flag)
%
% GRID_func()
% 
% DESCRIPTION:
% The function is a part of the GRID algortihm. This is part of the
% estimation section of the algortihm. The function finds the first
% derivative using the 3-point lagrangian interpolation polynomial. The
% first derivative is used later in the last section of the GRID algorithm.
% 
% INPUT:
% Delta G               - The maximum ROC (rate of change)            
% G                     - vector: dim: 3x1 ( [Gm-2, Gm-1, Gm] )
% tspan                 - The corresponding time values for the G vector
% prev_vec              - Vector of [G_{F,NS}(k-1), G_{F}(k-2)],eq: (1)&(3)
% Gmin                  - Vector of [G_{min,1},G_{min,2},G_{min,3}]
% Gfm_vec               - Vector of [G'_{F}(k-2),G'{F}(k-1)], eq: (4)
% flag                  - 
%
% OUTPUT:
% Gfm_vec               - The stored new vector of [G'_{F}(k-1),G'{F}(k-2)]           
% prev_vec              - The stored new vector of [G_{F,NS}(k-1), G_{F}(k-2)]
% zero_one              - 1 or 0 for detected meal.
% flag                  -
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

% Inisializing all values

Gfm_m2 = Gfm_vec(1);     % The second previous derivative used in euqation 4
Gfm_m1 = Gfm_vec(2);     % The previous derivative used in equation 4

Gfnsm2_prev = G_prev(1); % The second previous noise-spike filtered value
                         % used in equation 1   
Gfm2_prev   = G_prev(2); % The second previous low filterd value used in 
                         % equation 2

% Minimum values used in equation 4
Gmin1 = Gmin(1);         
Gmin2 = Gmin(2);         
Gmin3 = Gmin(3);            

% The two previous measured gluscose values and the one at control state
Gm2 = G(1);              % The second previous glucose value  
Gm1 = G(2);              % The previous glucose value  
G   = G(3);              % The glucose value at control state   


% COMPUTING 

% The noise-spike filter at the 3 sampling times
Gfnsm2 = spikefilt_func(Gm2,Gfnsm2_prev,delta_G);
Gfnsm1 = spikefilt_func(Gm1,Gfnsm2,delta_G);
Gfns   = spikefilt_func(G,Gfnsm1,delta_G);

% The low filter at the 3 sampling times
Gfm2 = lowfilt_func(tau,tspan,Gfnsm2,Gfm2_prev);
Gfm1 = lowfilt_func(tau,tspan,Gfnsm1,Gfm2);
Gf   = lowfilt_func(tau,tspan,Gfns,Gfm1);

% Inisializing input for estimate_lagrange
Gf_vec = [ Gfm2 , Gfm1 , Gf ]; 

% Computing the first derivative using lagrange 
Gfm    = estimate_lagrange(t_vec,Gf_vec); % Returns the derivative 

% The detection part from equation 4
if (Gf > Gmin1) && ... 
   ((Gfm > Gmin3) && (Gfm_m1 > Gmin3) && (Gfm_m2 > Gmin3) ...
   || (Gfm > Gmin2) && (Gfm_m2 > Gmin2) )
    
    zero_one = 1; % A meal has been detected 
    
else
    
    zero_one = 0; % No meal has been detected
    
end

% Output
G_prev = [Gfnsm2,Gfm2]; % Outputting the updated filtered values 
Gfm_vec = [Gfm_m1,Gfm]; % Outputting the updated derivative values 


% Counting part 
if flag > 0 
    
    flag = flag -1;
    zero_one = 0;

elseif flag == 0 && zero_one == 1
    
    flag = 120/tspan;
    
end 
   
    
end









