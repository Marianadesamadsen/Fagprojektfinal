   for i = 1:length(idx_missed)
    % Index value for missed bolus
    k_missed = idx_missed(i)
    
    % Index value for lessened bolus
    k_less = idx_less(i)
    
        if D(k_missed) == 1 && bolus_missed(k_missed) == 1 && sum(D_detected(k_missed:k_missed+stride)) == 0
            falsenegative = falsenegative + 1
        end
    
        if D(k_less) == 1 && bolus_less(k_less) == 1 && sum(D_detected(k_less:k_less+stride)) == 0
            falsenegative = falsenegative + 1
            
        end
   end
        