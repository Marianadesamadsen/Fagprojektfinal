

% 
% % %% Checking for false positve or false negative values
% 
% falsenegative   = 0;
% falsepositive   = 0;
% truepositive    = 0;
% count           = 0;
% 
% stride = 90 / Ts; 
% max    = 50 / Ts;
% 
% i = 1;
% 
% while i < length(zero_one) - stride
% 
%     if D(1,i) < max && sum(zero_one(i:i+stride)) == 1 && sum(D(1,i:i+stride)) < max 
%         falsepositive = falsepositive + 1; 
%         i = i + stride;
%         
%     elseif D(1,i) >= max && sum(zero_one(i:i+stride)) == 0 
%         falsenegative = falsenegative + 1;
%         i = i + stride;
%         
%     elseif D(1,i) >= max && sum(zero_one(i:i+stride)) == 1
%        truepositive = truepositive + 1;
%        i = i + stride;
%         
%     else
%         i = i + 1;
%     end 
%     
%     
% end
%     
% falsepositive1 = falsepositive;
% falsenegative1 = falsenegative;
% truepositive1 = truepositive;
% fprintf('number of false positive: %d \n',falsepositive1);
% fprintf('number of false negative: %d\n',falsenegative1);
% fprintf('number of true positive: %d\n',truepositive1);
