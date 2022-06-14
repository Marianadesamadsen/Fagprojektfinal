% test 
a=D(:,3)';
h=D_detected(:,1,3);

% Initializing 
falsenegative   = 0;
falsepositive   = 0;
truepositive    = 0;
truemeals       = zeros(1,N);
mealdetec       = zeros(1,N);

% Changing datatype of D to binary
for i = 1:N
    
    if a(1,i) >= 50/Ts % Not considering the snackmeals 
        a(1,i) = 1;
        
    else
        a(1,i) = 0;
    end 
    
end 

% Finding the indices for the truemeals
for i = 1 : N
    if a(1,i) == 1
        truemeals(i) = i;
    end 
end

% Finding the indices for the detected meals 
for i = 1 : N
    if h(i) == 1
        mealdetec(i) = i;
    end 
end

% Removing all the zeros so there is only the indices left
idxdetecmeals = nonzeros(mealdetec)';
idxtruemeals = nonzeros(truemeals)';

% Examine if there are no detected meals where there are true meals
for i = 1:length(idxtruemeals)
   
    % The idx value when there is a true meal
    k = idxtruemeals(i);
    
    if sum(h(k:k+stride)) == 0 
        falsenegative = falsenegative + 1;
    end  
end

% Examine if there are no true meals where there are detected meals 
for i = 1:length(idxdetecmeals)
    
    % The idx value when meal has been detected 
    k = idxdetecmeals(i);
    
    if (k-stride) < 1 
        j = k-1;
        
        if  sum(a(1,k-j:k)) == 0 
        falsepositive = falsepositive + 1;
        end
    
    elseif sum(a(1,k-stride:k)) == 0 
        falsepositive = falsepositive + 1;
    end 
    
end

% Examine if there are true meals where there are detected meals
for i = 1:length(idxdetecmeals)
   
    % The idx value when meal has been detected 
    k = idxdetecmeals(i);
    
    if (k-stride) < 1 
        j = k-1;
        
        if  sum(a(1,k-j:k)) == 0 
        falsepositive = falsepositive + 1;
        end
    
    elseif sum(a(1,k-stride:k)) == 1 
        truepositive = truepositive + 1;
    end  
end

% FINDING true negatives 

% Since TOTAL = FP + FN + TP + TN -> TN = TOTAL - FN - FP - TP
truenegative = N-falsenegative-falsepositive-truepositive;
