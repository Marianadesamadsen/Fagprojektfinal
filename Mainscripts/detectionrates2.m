function [truepositive, falsepositive, falsenegative, truenegative] = detectionrates2(stride,D,D_detected,Ts,idx_missed, idx_less,U)
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
%
% D                 - The true vector of real meals.
% D_detected        - The estimated vecor of 0 or 1. 1 meaning meals is
% detected
% U                 - Bolus insulin
%
% Ts                - The time between control steps
% idx_missed        - Indicies for missed bolus when ingesting a meal
% idx_less          - Indicies for lessened bolus when ingesting a meal
% U                 - Insulin matrix
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

% bolus vector
bolus = U(2,:);

% bolus missed and bolus lessened vectors such that we have a 8640 binary
% vector with entry 1 if there's missed or lessened bolus
bolus_missed = zeros(1,length(bolus));

    for i = 1 : length(idx_missed)
         bolus_missed(idx_missed(i)) = 1;
    end

bolus_less = zeros(1,length(bolus));

    for i = 1:length(idx_less)
        bolus_less(idx_less(i)) = 1;
    end

% Changing datatype of D to binary with entry 1 for a consumed meal
    for i = 1:length(D(1,:))

        if D(1,i) >= 50/Ts % Not considering the snackmeals
        D(1,i) = 1;

        else
        D(1,i) = 0;
        end

    end

% Changing datatype bolus to binary with entry 1 for given bolus
bolus(find(bolus)) = 1;

% Indicies for detected meal and true meal
idx_detected = find(D_detected);
idx_true_meal = find(D);

 % **** TRUE POSITIVE ****
    % If there's a detected true meal in the same stride range where either
    % bolus is missed or bolus is lessened

    for i = 1:length(idx_missed)

    % Index value for missed bolus
    k_missed = idx_missed(i);

    % Index value for lessened bolus
    k_less = idx_less(i);

    % Index value for D true meal

    % If there's detected a true meal where bolus is missed
    if D(k_missed)==1 && sum(D_detected(k_missed:k_missed+stride)) == 1
        truepositive = truepositive + 1;
    end

    % If there's a detected a true meal where bolus is lessened
    if D(k_less) == 1 && sum(D_detected(k_less:k_less+stride)) == 1
        truepositive = truepositive + 1;
    end

    end

    % **** FALSE POSITIVE ****
    % If there's detected meals when there's been given bolus insulin

    for i = 1:length(nonzeros(D_detected))
       % indicies for detected meals
       k_detected = idx_detected(i);

       % Considering both true detected meal and neither missed nor
       % lessened bolus
       if D_detected(k_detected) == 1 && sum(D(k_detected-stride:k_detected)) == 1 && sum(bolus_missed(k_detected-stride:k_detected)) == 0 && sum(bolus_less(k_detected-stride:k_detected)) == 0
           falsepositive = falsepositive + 1;
       end
       if D_detected(k_detected) == 1 && sum(D(k_detected-stride:k_detected)) == 0
           falsepositive = falsepositive + 1;
       end

    end

    % **** FALSE NEGATIVE ****
    % When bolus is either missed or lessened but no meal is detected

  for i = 1:length(idx_missed)
    % Index value for missed bolus
    k_missed = idx_missed(i);

    % Index value for lessened bolus
    k_less = idx_less(i);

        if D(k_missed) == 1 && bolus_missed(k_missed) == 1 && sum(D_detected(k_missed:k_missed+stride)) == 0
            falsenegative = falsenegative + 1;
        end

        if D(k_less) == 1 && bolus_less(k_less) == 1 && sum(D_detected(k_less:k_less+stride)) == 0
            falsenegative = falsenegative + 1;
        end
  end
   
  truenegative = length(D) - truepositive - falsepositive - falsenegative;


end
