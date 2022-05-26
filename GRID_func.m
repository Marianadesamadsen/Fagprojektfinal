function  [ Gfm_vec , G_prev , flag, zero_one ] = ...
            GRID_func( delta_G , G , tau, tspan , G_prev , Gmin, Gfm_vec , t_vec, flag)
%        
% GOAL:
% Making the GRID function, which outputs 0 or 1 depending on having
% detected a meal.
%
% INPUT:
% Delta G:              The maximum ROC (rate of change)            
% G:                    vector: dim: 3x1 ( [Gm-2, Gm-1, Gm] )
% tspan:                The corresponding time values for the G vector
% prev_vec:             Vector of [G_{F,NS}(k-1), G_{F}(k-2)],eq: (1)&(3)
% Gmin:                 Vector of [G_{min,1},G_{min,2},G_{min,3}]
% Gfm_vec:              Vector of [G'_{F}(k-2),G'{F}(k-1)], eq: (4)
% 
% OUTPUT:   
% Gfm_vec:              The stored new vector of [G'_{F}(k-1),G'{F}(k-2)]           
% prev_vec:             The stored new vector of [G_{F,NS}(k-1), G_{F}(k-2)]
% zero_one:             1 or 0 for detected meal.
% 

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



% Inizialising
Gfm_m2 = Gfm_vec(1);
Gfm_m1 = Gfm_vec(2);

Gfnsm2_prev = G_prev(1); % eq 1
Gfm2_prev   = G_prev(2); % eq 2

Gmin1 = Gmin(1);
Gmin2 = Gmin(2);
Gmin3 = Gmin(3);

Gm2 = G(1); 
Gm1 = G(2); 
G   = G(3); 

% First filter function
Gfnsm2 = spikefilt_func(Gm2,Gfnsm2_prev,delta_G);
Gfnsm1 = spikefilt_func(Gm1,Gfnsm2,delta_G);
Gfns   = spikefilt_func(G,Gfnsm1,delta_G);

% Second filter function
Gfm2 = lowfilt_func(tau,tspan,Gfnsm2,Gfm2_prev);
Gfm1 = lowfilt_func(tau,tspan,Gfnsm1,Gfm2);
Gf   = lowfilt_func(tau,tspan,Gfns,Gfm1);

% Lagrange
Gf_vec = [ Gfm2 , Gfm1 , Gf ]; 

Gfm    = estimate_lagrange(t_vec,Gf_vec);

% Grid 
if (Gf > Gmin1) && ( (Gfm > Gmin3) && (Gfm_m1 > Gmin3) ...
        && (Gfm_m2 > Gmin3) || (Gfm > Gmin2) && (Gfm_m2 > Gmin2) )
    
    zero_one = 1; % detected meal
    
else
    
    zero_one = 0; % no meal detected
    
end

% Output
G_prev = [Gfnsm2,Gfm2]; % The previous for: Gfns-1 and Gf-1 eq. 1 & 2
Gfm_vec = [Gfm_m1,Gfm]; % The previous for: Gfm-m1 and Gfm-m2

if flag > 0 
    
    flag = flag -1;
    zero_one = 0;

elseif flag == 0 && zero_one == 1
    
    flag = 120/tspan;
    
end 
   
    
end









