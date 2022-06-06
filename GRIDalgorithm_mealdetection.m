function D_detected = GRIDalgorithm_mealdetection(G,Gmin,tau,delta_G,t_vec,Ts)

% Inisializing 
filt_prev      = zeros(length(G),2); % The vector of previous filtered values
Gfm_vec        = zeros(length(G),2); % The vector of previous derivatives
G_vec          = [G(1),G(1),G(1)];   % Inserting the start previous glucose measurements as the same value
filt_prev(1,:) = [G(1),G(1)];        % Inserting the previous filtered value as the not filtered values
flag           = 0;                  % No detected meal to begin with


% Computing the first two detections
[ Gfm_vec(2,:) , filt_prev(2,:) , flag, D_detected(2) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(1,:) , Gmin, Gfm_vec(1,:) , t_vec ,flag );
                                
                                % Updating G_vec
                                G_vec=[G(1),G(1),G(2)];
                                
[ Gfm_vec(3,:) , filt_prev(3,:) , flag, D_detected(3) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(2,:) , Gmin, Gfm_vec(2,:) , t_vec, flag );
                                
                                % Updating G_vec
                                G_vec=[G(1),G(2),G(3)];
                                
% Computing the last detections
for i = 3 : length(G)-1
    
[ Gfm_vec(i+1,:) , filt_prev(i+1,:) , flag, D_detected(i) ] = GRID_func( delta_G , G_vec , tau, Ts , ...
                                    filt_prev(i,:) , Gmin, Gfm_vec(i,:) , t_vec , flag );
                                
                                % Updating G_vec
                                G_vec=[G(i-1),G(i),G(i+1)];
                                
end

end

