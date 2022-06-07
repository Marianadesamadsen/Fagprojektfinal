function [truenegative,truepositive,falsepositive,falsenegative] = detectionrates(stride,D,D_detected,Ts)

falsenegative   = 0;
falsepositive   = 0;
truepositive    = 0;
truenegative    = 0;
truemeals       = zeros(1,length(D(1,:)));
mealdetec       = zeros(1,length(D_detected));

% Chaning datatype to binary
for i = 1:length(D(1,:))
    
    if D(1,i) >= 50/Ts
        D(1,i) = 1;
        
    else
        D(1,i) = 0;
    end 
    
end 

% Finding the indices for the truemeals
for i = 1 : length(D(1,:))
    if D(1,i) == 1
        truemeals(i) = i;
    end 
end

% Finding the indices for the detected meals
for i = 1 : length(D(1,:))
    if D_detected(i) == 1
        mealdetec(i) = i;
    end 
end

% Removing all the zeros
idxdetecmeals = find(mealdetec);
idxtruemeals = find(truemeals);

% Examine if there are no detected meals where there are true meals
for i = 1:length(idxtruemeals)
   
    k = idxtruemeals(i);
    
    if sum(D_detected(k:k+stride)) == 0 
        falsenegative = falsenegative + 1;
    end  
end

% Examine if there are no true meals where there are detected meals 
for i = 1:length(idxdetecmeals)
   
    k = idxdetecmeals(i);
    
    if sum(D(1,k-stride:k)) == 0 
        falsepositive = falsepositive + 1;
    end  
end

% Examine if there are true meals where there are detected meals
for i = 1:length(idxdetecmeals)
   
    k = idxdetecmeals(i);
    
    if sum(D(1,k-stride:k)) == 1 
        truepositive = truepositive + 1;
    end  
end

% FINDING true negatives 
% 
% for i = 1:length(D(1,:))
%     
%     if D(1,i) == 0 && D_detected(i) == 0 
%        truenegative = truenegative + 1;
%     end 
%     
% end 
% 
% % Adding the truepositive since the code does not consider that the
% % detected meals are staggered from the actual meals. 
% truenegative = truenegative + truepositive;

truenegative = length(D(1,:))-falsenegative-falsepositive-truepositive;

end


