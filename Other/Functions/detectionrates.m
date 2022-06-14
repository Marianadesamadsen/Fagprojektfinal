function [truenegative,truepositive,falsepositive,falsenegative] = detectionrates(stride,D,D_detected,Ts)
%
% detectionrates()
% 
% DESCRIPTION:
% This function computes the TP, FP, TN, FN based on a true vector with
% detected meals D, and a computed vector with detected meals D_detected. 
% It uses a stride since the meal will be detected a time period after the 
% meals was given.
% It uses the indices for the true meals and the detected meals to compare
% if in the stride a meals should be detected or not and vise versa. 
%
% INPUT:
% stride            - The maximal time period it takes from the meal to be
% detected 
% D                 - The true vector of real meals.
% D_detected        - The estimated vecor of 0 or 1. 1 meaning meals is
% detected
% U                 - Bolus insulin
% 
% Ts                - The time between control steps
%
% OUTPUT:
% Four outputs being TP, FP, TN, FN
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

% Initializing 
falsenegative   = 0;
falsepositive   = 0;
truepositive    = 0;
truemeals       = zeros(1,length(D(1,:)));
mealdetec       = zeros(1,length(D_detected));

% Changing datatype of D to binary
for i = 1:length(D(1,:))
    
    if D(1,i) >= 50/Ts % Not considering the snackmeals 
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

% Removing all the zeros so there is only the indices left
idxdetecmeals = nonzeros(mealdetec)';
idxtruemeals = nonzeros(truemeals)';

% Examine if there are no detected meals where there are true meals
for i = 1:length(idxtruemeals)
   
    % The idx value when there is a true meal
    k = idxtruemeals(i);
    
    if sum(D_detected(k:k+stride)) == 0 
        falsenegative = falsenegative + 1;
    end  
end

% Examine if there are no true meals where there are detected meals 
for i = 1:length(idxdetecmeals)
    
    % The idx value when meal has been detected 
    k = idxdetecmeals(i);
    
    if sum(D(1,k-stride:k)) == 0 
        falsepositive = falsepositive + 1;
    end  
end

% Examine if there are true meals where there are detected meals
for i = 1:length(idxdetecmeals)
   
    % The idx value when meal has been detected 
    k = idxdetecmeals(i);
    
    if sum(D(1,k-stride:k)) == 1 
        truepositive = truepositive + 1;
    end  
end

% FINDING true negatives 

% Since TOTAL = FP + FN + TP + TN -> TN = TOTAL - FN - FP - TP
truenegative = length(D(1,:))-falsenegative-falsepositive-truepositive;

end


