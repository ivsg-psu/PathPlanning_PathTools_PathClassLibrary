function fcn_medianFilterViaSplineFitting2(xEast,xEast_increment,xEast_increment_pred,DGPS_is_active)
         

%fcn_medianFilterViaSplineFitting2 - smooths the xEast values.
% Step 1: uses merged and smooted values for increments (from all sensors)
% to calculate where the jump locations are in the current sensor. 

% Revision history:
% 2019_11_27 - first write of function, moving material out of main code
% area.

threshold = 10;
flag_plot_variances = 1;

% Find the locations where DGPS is on
index_DGPS_is_active = find(DGPS_is_active == 1);

% Find the locations where DGPS is on at the start
index_DGPS_shuts_off_first = find(diff(DGPS_is_active)==-1,1) -1;
start_good_indices = 1:index_DGPS_shuts_off_first;

% Find the locations where DGPS is on at the end
index_DGPS_turns_on_last = find(diff(DGPS_is_active)==1,1,'last') + 1;
end_good_indices = index_DGPS_turns_on_last:length(DGPS_is_active(:,1));

% Find the locations of the middle indices
middle_bad_indices = (index_DGPS_shuts_off_first+1):(index_DGPS_turns_on_last-1);

% Find the locations where jumps in bias occur, and cumulative jumps
[index_jump_locations,cum_jumps] = fcn_calculateJumpLocations(xEast,xEast_increment_pred);

% Calculate differences between predicted and actual incrments at start
mean_bias_start = mean(cum_jumps(start_good_indices,1));
cum_jumps = cum_jumps - mean_bias_start; % Shift cumulative jumps by this mean bias, to force them to start to zero
mean_bias_start = mean(cum_jumps(start_good_indices,1)); % Recalculate differences between predicted and actual incrments at start

% Now see what the bias is at the end...
mean_bias_end = mean(cum_jumps(end_good_indices,1));
bias_drift = mean_bias_end - mean_bias_start;
bias_change = cum_jumps(index_DGPS_turns_on_last) - cum_jumps(index_DGPS_shuts_off_first);

% Calculate the scaling at the middle to remove the bias_drift

fprintf(1,'BEFORE FIXING:');
fprintf(1,'\tAve start bias: %f meters\n',mean_bias_start);
fprintf(1,'\tBias when DGPS shut off (should be close to zero): %f: meters\n',cum_jumps(index_DGPS_shuts_off_first))
fprintf(1,'\tAve end bias: %f meters\n',mean_bias_end);
fprintf(1,'\tBias when DGPS turns on at end: %f: meters\n',cum_jumps(index_DGPS_turns_on_last))
fprintf(1,'\tDrift in bias from averages: %f meters\n',bias_drift);
fprintf(1,'\tChange in bias calculated during outage: %f meters\n',bias_change);

% Initialize the new bias vector
new_cum_jumps = cum_jumps;

% Calculate offsets for first jumps such that it cancels the start bias
goal_trend_indices = [start_good_indices(1); start_good_indices(end)];
goal_trend_values  = [0; -mean_bias_start];
offset_to_cancel_bias = interp1(goal_trend_indices,goal_trend_values,start_good_indices);
new_cum_jumps(start_good_indices) = new_cum_jumps(start_good_indices) + offset_to_cancel_bias';

% Calculate offsets for middle jumps such that it cancels the end bias
goal_trend_indices = [middle_bad_indices(1); middle_bad_indices(end)];
goal_trend_values  = [0; -mean_bias_end];
% Format is: Vq = interp1(X,V,Xq
offset_to_cancel_bias = interp1(goal_trend_indices,goal_trend_values,middle_bad_indices);
new_cum_jumps(middle_bad_indices) = new_cum_jumps(middle_bad_indices) + offset_to_cancel_bias';

% Calculate offsets for the end area such that it cancels the end bias
offset_to_cancel_bias = -mean_bias_end*ones(1,length(end_good_indices));
new_cum_jumps(end_good_indices) = new_cum_jumps(end_good_indices) + offset_to_cancel_bias';

figure(28282); clf; hold on; grid minor; plot(offset_to_cancel_bias)

% RECALCULATE differences between predicted and actual incrments at start
mean_bias_start = mean(new_cum_jumps(start_good_indices,1));
new_cum_jumps = new_cum_jumps - mean_bias_start; % Shift cumulative jumps by this mean bias, to force them to start to zero
mean_bias_start = mean(new_cum_jumps(start_good_indices,1)); % Recalculate differences between predicted and actual incrments at start

% RECALCULATE to see what the bias is at the end...
mean_bias_end = mean(new_cum_jumps(end_good_indices,1));
bias_drift = mean_bias_end - mean_bias_start;
bias_change = new_cum_jumps(index_DGPS_turns_on_last) - new_cum_jumps(index_DGPS_shuts_off_first);

% Calculate the scaling at the middle to remove the bias_drift

fprintf(1,'AFTER FIXING:');
fprintf(1,'\tAve start bias: %f meters\n',mean_bias_start);
fprintf(1,'\tBias when DGPS shut off (should be close to zero): %f: meters\n',new_cum_jumps(index_DGPS_shuts_off_first))
fprintf(1,'\tAve end bias: %f meters\n',mean_bias_end);
fprintf(1,'\tBias when DGPS turns on at end: %f: meters\n',new_cum_jumps(index_DGPS_turns_on_last))
fprintf(1,'\tDrift in bias from averages: %f meters\n',bias_drift);
fprintf(1,'\tChange in bias calculated during outage: %f meters\n',bias_change);





xEast_clean_no_jumps = xEast - new_cum_jumps;

figure(3444);
clf
hold on;
plot(xEast,xEast_increment,'k');
plot(xEast,xEast_increment_pred,'m');
plot(xEast(index_DGPS_is_active),xEast_increment(index_DGPS_is_active),'b.','LineWidth',3); %#ok<*FNDSB>
plot(xEast(index_jump_locations),xEast_increment(index_jump_locations),'ro');
plot(xEast_clean_no_jumps,xEast_increment_pred,'g');


return

function [index_jump_locations,jump_mag] = fcn_calculateJumpLocations(x,x_increment_predicted)
flag_make_plots = 0;
jump_limit = 0.1;
N_standard_deviations = 3;

% calculate the biases and change in biases
x_increments = [0; diff(x(:,1))];
biases = cumsum(x_increments) - cumsum(x_increment_predicted);
change_in_biases = [0;diff(biases)];

% The standard deviation will be grossly skewed by jumps, so calculate a
% "clean" standard deviationt hat does not have the jumps within the
% population. 
indices_good_change_in_biases = find(change_in_biases<jump_limit & change_in_biases > -jump_limit);
std_biases = std(change_in_biases(indices_good_change_in_biases)); 

% For debugging:
% figure(3344);clf;histogram(change_in_biases,10000);
% title(sprintf('Standard deviation in biases, with outliers removed, is: %f',std_biases));

% Try to calculate the jumps by finding locations where they jump larger
% than N number of standard deviations.
index_jump_locations = find(abs(change_in_biases)>N_standard_deviations*std_biases);

% Fill in the jumps
jumps = 0*biases;  % Initialize the vector
jumps(index_jump_locations) = change_in_biases(index_jump_locations);
jump_mag = cumsum(jumps);
true_bias_drift = biases - jump_mag;

if 1==flag_make_plots
    figure(2);
    clf;
    hold on;
    time = (1:length(biases))'*0.05;
    plot(time, biases);
    plot(time(index_jump_locations),biases(index_jump_locations),'ro');
    plot(time,true_bias_drift,'m');
end
return

