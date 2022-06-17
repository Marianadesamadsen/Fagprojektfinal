function D_detected = GRIDalgorithm_mealdetection2(G,Gmin,tau,delta_G,tspan,Ts)
% GRIDalgorithm_mealdetection()
% 
% DESCRIPTION:
% This function takes in glucose measurements (G), that consists of several
% control steps and computes the whole meal detection part of the GRID
% algorithm. This means it is extended to several datapoints not only 1. It
% differs from GRIDalgorithm_mealdetection by being able to test the
% GRID-algorithm on any given glucose measurements where the control steps
% are not equally distanced. 
%
% INPUT:
% G             - The glucose concentration for several datapoints
% Gmin          - Vector of Gmin values for the detection part
% tau           - Filter time constant 
% delta_ G      - Maximum allowable ROC (from article)
% t_vec         - The current time and the two previous time points 
% Ts            - The step size between control steps. 
%
% OUTPUT:
% D_detected    - The vector with 0 or 1. 1 meaning meal is detected. 
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


% Inisializing 
filt_prev      = zeros(length(G),2); % The vector of previous filtered values
Gfm_vec        = zeros(length(G),2); % The vector of previous derivatives
G_vec          = [G(1),G(1),G(1)];   % Inserting the start previous glucose measurements as the same value
filt_prev(1,:) = [G(1),G(1)];        % Inserting the previous filtered value as the not filtered values
flag           = 0;                  % No detected meal to begin with
t_vec          = [tspan(1),tspan(1),tspan(1)];
Ts_t           = Ts(1);



% Computing the first two detections
[ Gfm_vec(2,:) , filt_prev(2,:) , flag, D_detected(2) ] = GRID_func( delta_G , G_vec , tau, Ts_t , ...
                                    filt_prev(1,:) , Gmin, Gfm_vec(1,:) , t_vec ,flag );
                                
                                % Updating G_vec and t_vec and Ts
                                G_vec=[G(1),G(1),G(2)];
                                t_vec=[tspan(1),tspan(1),tspan(2)];
                                Ts_t=Ts(2);
                               
[ Gfm_vec(3,:) , filt_prev(3,:) , flag, D_detected(3) ] = GRID_func( delta_G , G_vec , tau, Ts_t , ...
                                    filt_prev(2,:) , Gmin, Gfm_vec(2,:) , t_vec, flag );
                                
                                % Updating G_vec and t_vec and Ts
                                G_vec=[G(1),G(2),G(3)];
                                t_vec=[tspan(1),tspan(2),tspan(3)];
                                Ts_t=Ts(3);
                                
% Computing the last detections
for i = 3 : length(G)-1
    
[ Gfm_vec(i+1,:) , filt_prev(i+1,:) , flag, D_detected(i) ] = GRID_func( delta_G , G_vec , tau, Ts_t , ...
                                    filt_prev(i,:) , Gmin, Gfm_vec(i,:) , t_vec , flag );
                                
                                % Updating G_vec and t_vec
                                G_vec=[G(i-1),G(i),G(i+1)];
                                t_vec=[tspan(i-1),tspan(i),tspan(i+1)];
                                Ts_t=Ts(i+1);
                                
end

end