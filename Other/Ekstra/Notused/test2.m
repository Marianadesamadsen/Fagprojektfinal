for i = 1:length(idx_missed)
    
    % Index value for missed bolus
    k_missed = idx_missed(i)
    
    % Index value for lessened bolus
    k_less = idx_less(i)
    
    % **** TRUE POSITIVE ****
    % If there's a detected true meal in the stride-range around the missed
    % bolus, the true positive is one
    if sum(D_detected(k_missed:k_missed+stride)) == 1 && sum(D(k_missed:k_missed+stride)) == 1
        truepositive = truepositive + 1
    end
    
    % If there's a detected true meal in the stride-range around the
    % lessened bolus, the true positive is one
    if sum(D_detected(k_less:k_less+stride)) == 1 && sum(D(k_missed:k_missed+stride)) == 1
        truepositive = truepositive + 1
    end
end
    
  for i = 1:length(D)
        k_true = idx_true_meal(i)
        if sum(D_detected(k_true-stride:k_true+stride)) == 0 && sum(bolus_missed(k_true-stride:k_true+stride))==0 && sum(bolus_less(k_true-stride:k_true+stride))==0
        truenegative = truenegative + 1
        end
    end