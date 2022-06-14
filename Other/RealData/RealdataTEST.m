% %% Test for data
%
% date_test = (t2(70:85));
% for i = 1:length(date_test)
%     if isequal(date_test(i), date_test(i+1)) == 1
%         date_test(i) = []
%     end
% end

% %% ***** FUNGERER IKKE The non duplicated datetimes for insulin, no repeats of time****
% time_CGM = unique(t);
% time_insulin = unique(t2);
%
% % The most measurements within almost the same time interval we have is
% % from the CGM measurements
% time_range = length(time_CGM);
% % Adding 100 extra datetimes in each end of the time interval
% extendFactor = 20;
% % Insulin starts getting measured 2 minutes before CGM and finishes 1
% % minute after CGM
% t_1 = time_insulin(1);
% t_2 = time_insulin(end);
%
% time_vec = linspace(t_1, t_2, time_range *extendFactor);
%
% % Adjusting the bolus insulin vector to match the time_vec
% %vec_ratio = length(time_vec)/length(ubo);  % length ratio between time vector and bolus vector. Is gonna be used for indexing
% %zero_ubo_vec = zeros(1,length(time_vec));  % allocating space for the extended bolus vector
%
% idx_bolus = find(ubo);
%
% zero_ubo_vec(1,length(time_vec));
%
% for i = 1:length(idx_bolus)
%     zero_ubo_vec(idx_bolus(i)*extendFactor) = ubo(idx_bolus(i));
% end
%
% ubo_vec = zero_ubo_vec;

%% The non duplicated datetimes for insulin, no repeats of time
%time_CGM = unique(t);       % Filtered time for CGM
%time_insulin = unique(t2);  % Filtered time for insulin

%t_start = min(time_CGM(1),time_insulin(1)); % Start time, the one out of the two measurements that start first
%t_end = max(time_CGM(end), time_insulin(end)); % Finish time, the one out of the two measurements that finishes last
%
% % timerange = max(length(time_CGM),length(time_insulin));
%
% time_vec2 = [t_start, time_CGM, t_end];