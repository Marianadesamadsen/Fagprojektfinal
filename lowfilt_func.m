function Gf_ctrlstate = lowfilt_func(tau,tspan,Gfns,Gf_prev)
% FAGPROJEKT 2022
% 
% ARTHUER:
% Emma Lind, Mona Saleem, Mariana de Sá Madsen
%
% GOAL: 
% To filter high frequency fluctations
% 
% INPUT:
% tau:          Filter time constant
% delta_t:      The samling time period 
% Gf:           The filtered value from spikefilt_func
% Gf_prev:      The previous filtered balue from lowfilt_func
%
% OUTPUT:
% The low filtered value at the glucose state
% 

% The equation is from the GRID low filter function
Gf_ctrlstate = (tspan / ( tau + tspan)) * Gfns + ( 1 - (tspan / (tau + tspan) ) )...
       * Gf_prev ; 

end

