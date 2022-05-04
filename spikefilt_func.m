function gfns_ctrlstate = spikefilt_func(Gm,Gfns_prev,delta_G)
%
% GOAL: 
% To filter the extrem measurements 
%
% INPUT: 
% GM:       The measuret glucose in the blood 
% Gf_prev   The filteret measuret glucose in the blood rigt before
% delta_G   The maximum ROC (rate of change)
% 
% OUTPUT:
% The filtered value at the glucose state 

% 3 if-statements followed by the GRID algorithm (preprossing section)

    if abs( Gm - Gfns_prev ) <= delta_G
        
        % If the change in glucose is smaller than the maximum 
        % allowed change then it is the same
        
        gfns_ctrlstate = Gm; 
        
    elseif  (Gfns_prev - Gm) > delta_G
        
        % If the change is larger, then we subtract the maximum ROC
        
        gfns_ctrlstate = Gfns_prev - delta_G;
    
    elseif  (Gm - Gfns_prev) > delta_G
        
        % If the change is larger, then we subtract the maximum ROC
        
        gfns_ctrlstate = Gfns_prev + delta_G;
    
    end
   
    
end








