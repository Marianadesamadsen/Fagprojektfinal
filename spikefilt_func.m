function Gfns_ctrlstate = spikefilt_func(G,Gfns_prev,delta_G)
% 
% spikefilt_func()
% 
% DESCRIPTION:
% The function is a part of the GRID algortihm. This is part of the
% preprocessing section of the algortihm. The function filters the data
% using a noise-spike filter. This means the function will consider what
% data are extreme and what is not. Such that only extreme measurements of
% glucose will be consideret as meals in the GRID algorithm.
%
% INPUT:
% G          - The glucose measurement in the blood at the control state
% Gfns_prev  - The previous filteret measuret glucose in the blood
% delta_G    - The maximum ROC (rate of change)
%
% OUTPUT:
% The filtered value of glucose measured at the control state
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


% 3 if-statements in the preprocessing section to detect extreme
% measurements

    if abs( G - Gfns_prev ) <= delta_G
        
        % If the change in glucose is smaller than the maximum 
        % rate of change then the measurement is not changed
        
        Gfns_ctrlstate = G; 
        
    elseif  (Gfns_prev - G) > delta_G
        
        % If the above change is larger, then maximum ROC is subtracted 
        % from the measurement
        
        Gfns_ctrlstate = Gfns_prev - delta_G;
    
    elseif  (G - Gfns_prev) > delta_G
        
        % If the above change is larger, then maximum ROC is added 
        % to the measurement
        
        Gfns_ctrlstate = Gfns_prev + delta_G;
    
    end
   
    
end








