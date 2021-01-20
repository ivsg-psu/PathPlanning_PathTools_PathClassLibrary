function [x_clean,x_increments_clean] = fcn_medianFilterViaIncrementTemplate(x,x_increment_pred,DGPS_is_active)
%fcn_medianFilterViaIncrementTemplate - smooths the state values by
% using a predicted trajectory in state space, provided by an
% increment_prediction, to that of the state, matching the average offsets
% in the locations at start and end of the trajectory where DGPS is active.
% 
% Step 1: using the data in vector DGPS_is_active, mark the data where DGPS
% corrections should be trusted. Split these into 3 sections: start,
% middle, end, where the start and end are good data, e.g. locked in, but
% middle may have drift.
% 
% Step 2: compare the predicted and actual trajectory, finding locations
% where the differences jump by a large number of standard deviations (4).
% At these jump locations, accumulate the jumps as a cumulative set of
% biases specifically caused by jumps.
%
% Step 3: Calculate how jump biases, if applied, skew data at start and
% end where DGPS should be trusted. 
%
% Step 4: remove the skew throughout start, middle, and end of trajectory
% by assuming a constant skew bias effect throughout start, middle, or end.
% Add this skew effect to the jump biases.
%
% Step 5: apply jump bias offsets to the state vector, which should now be
% smooth even in cases where no DGPS, and with near zero error in matching
% start/end of DGPS-active regions.

% Revision history:
% 2019_11_27 - first write of function, moving material out of main code
% area.
% 2019_12_01 - finished editing. Functionalized the code.

flag_plot_results = 0;
flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\n\t\tWithin subfunction: %s\n',namestr);
    fprintf(1,'\t\t\tUsing templating and median filtering to correct differntial jumps.\n');   
    fprintf(1,'\t\t\tLength of X: %d\n',length(x(:,1)));
    fprintf(1,'\t\t\tLength of predicted X increments: %d\n',length(x_increment_pred(:,1)));  
end

%% Check input data
if any(isnan(x)) || any(isnan(x_increment_pred)) || any(isnan(DGPS_is_active))
    error('NaN values detected on input data - unable to continue.');
end

if length(x(:,1))~=length(x_increment_pred(:,1))
   
    error('X and predicted X_increments need to have the same length');
end

if any(DGPS_is_active(1:50,1)==0)
    error('DGPS data must exist at start of this time slice for at least 50 samples');
end

if any(DGPS_is_active(end-50:end,1)==0)
    error('DGPS data must exist at end of this time slice for at least 50 samples');
end

%% Generate the increment vector (useful later)
x_increment = diff(x);
x_increment = [x_increment_pred(1,1); x_increment];

% Shift the predicted increments to make sure means match 
% (note: this does not affect final results at all, as these are based on
% differences only)
offset = mean(x_increment_pred) - mean(x_increment);
x_increment_pred = x_increment_pred - offset;

%% Step 1: mark start, middle, and end locations of data based on DGPS lock 
% Note: we add a buffer to each of at least 1 second, or 20 samples
num_buffer_indices  = 50;

% Find the locations where DGPS is on at the start
index_DGPS_shuts_off_first = find(diff(DGPS_is_active)==-1,1) -num_buffer_indices;
start_good_indices = 1:index_DGPS_shuts_off_first;

% Find the locations where DGPS is on at the end
index_DGPS_turns_on_last = find(diff(DGPS_is_active)==1,1,'last') + num_buffer_indices;
end_good_indices = index_DGPS_turns_on_last:length(DGPS_is_active(:,1));

% Find the locations of the middle indices
middle_bad_indices = (index_DGPS_shuts_off_first+1):(index_DGPS_turns_on_last-1);

%% Step 2: Calculate where jumps occur in prediction versus measurement, and find cumulative bias
% Find the locations where jumps in bias occur, and cumulative jumps
[index_jump_locations,cum_jumps] = fcn_calculateJumpLocations(x,x_increment_pred);

%% Step 3: Calculate how jump biases, if applied, skew data at start and end
% Calculate differences between predicted and actual increments at start
mean_bias_start = mean(cum_jumps(start_good_indices,1));
cum_jumps = cum_jumps - mean_bias_start; % Shift cumulative jumps by this mean bias, to force them to start to zero
mean_bias_start = mean(cum_jumps(start_good_indices,1)); % Recalculate differences between predicted and actual incrments at start

% Now see what the bias is at the end...
mean_bias_end = mean(cum_jumps(end_good_indices,1));
bias_drift = mean_bias_end - mean_bias_start;
bias_change = cum_jumps(index_DGPS_turns_on_last) - cum_jumps(index_DGPS_shuts_off_first);

% Calculate the scaling at the middle to remove the bias_drift
if 1==flag_do_debug
    fprintf(1,'\t\t\tBEFORE FIXING:');
    fprintf(1,'\t\t\t\tAve start bias: %f meters\n',mean_bias_start);
    fprintf(1,'\t\t\t\tBias when DGPS shut off (should be close to zero): %f: meters\n',cum_jumps(index_DGPS_shuts_off_first))
    fprintf(1,'\t\t\t\tAve end bias: %f meters\n',mean_bias_end);
    fprintf(1,'\t\t\t\tBias when DGPS turns on at end: %f: meters\n',cum_jumps(index_DGPS_turns_on_last))
    fprintf(1,'\t\t\t\tDrift in bias from averages: %f meters\n',bias_drift);
    fprintf(1,'\t\t\t\tChange in bias calculated during outage: %f meters\n',bias_change);
end

%% Step 4: remove the skew throughout start, middle, and end of trajectory
% by assuming a constant skew bias effect throughout start, middle, or end.
% Add this skew effect to the jump biases.

% Initialize the new bias vector
new_cum_jumps = cum_jumps;

flag_remove_skew = 1;
if flag_remove_skew
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
end
% For debugging
%figure(28282); clf; hold on; grid minor; plot(offset_to_cancel_bias)

% RECALCULATE differences between predicted and actual incrments at start
mean_bias_start = mean(new_cum_jumps(start_good_indices,1));
new_cum_jumps = new_cum_jumps - mean_bias_start; % Shift cumulative jumps by this mean bias, to force them to start to zero
mean_bias_start = mean(new_cum_jumps(start_good_indices,1)); % Recalculate differences between predicted and actual incrments at start

% RECALCULATE to see what the bias is at the end...
mean_bias_end = mean(new_cum_jumps(end_good_indices,1));
bias_drift = mean_bias_end - mean_bias_start;
bias_change = new_cum_jumps(index_DGPS_turns_on_last) - new_cum_jumps(index_DGPS_shuts_off_first);

% Calculate the scaling at the middle to remove the bias_drift
if 1==flag_do_debug
    
    fprintf(1,'\t\t\tAFTER FIXING:');
    fprintf(1,'\t\t\t\tAve start bias: %f meters\n',mean_bias_start);
    fprintf(1,'\t\t\t\tBias when DGPS shut off (should be close to zero): %f: meters\n',new_cum_jumps(index_DGPS_shuts_off_first))
    fprintf(1,'\t\t\t\tAve end bias: %f meters\n',mean_bias_end);
    fprintf(1,'\t\t\t\tBias when DGPS turns on at end: %f: meters\n',new_cum_jumps(index_DGPS_turns_on_last))
    fprintf(1,'\t\t\t\tDrift in bias from averages: %f meters\n',bias_drift);
    fprintf(1,'\t\t\t\tChange in bias calculated during outage: %f meters\n',bias_change);
end



%% Step 5: apply jump bias offsets to the state vector, which should now be
% smooth even in cases where no DGPS, and with near zero error in matching
% start/end of DGPS-active region

x_clean = x - new_cum_jumps;
x_increments_clean = x_increment;
x_increments_clean(index_jump_locations) = x_increment_pred(index_jump_locations);

if flag_plot_results
    figure(3424); clf;  hold on; grid minor;
    plot(x,x_increment,'k');
    xlabel('x (m)');
    ylabel('x increments (m)');
    title('Original data');
    legend('Raw data');
    
    figure(3524); clf;  hold on; grid minor;
    plot(x,x_increment,'k');
    plot(x(index_jump_locations),x_increment(index_jump_locations),'ro');
    xlabel('x (m)');
    ylabel('x increments (m)');
    title('Using the jumps to find discontinuities');
    legend('Raw data','Locations of jumps');
    
    figure(3624); clf;  hold on; grid minor;
    plot(x,x_increment,'k');
    plot(x(index_jump_locations),x_increment(index_jump_locations),'ro');
    plot(x(1:index_DGPS_shuts_off_first,1),x_increment(1:index_DGPS_shuts_off_first,1),'bd','MarkerSize',4,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0.6 0.6 1]); 
    plot(x(index_DGPS_turns_on_last:end,1),x_increment(index_DGPS_turns_on_last:end,1),'bd','MarkerSize',4,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0.6 0.6 1]);
    plot(x,x_increment_pred,'m');
    xlabel('x (m)');
    ylabel('x increments (m)');
    title('Comparing predicted template to original to find skew.');
    legend('Raw data','Locations of jumps','DGPS lock at start','DGPS lock at end','Predicted data');
         
    figure(3724); clf;  hold on; grid minor;
    plot(x,x_increment,'k');
    plot(x(index_jump_locations),x_increment(index_jump_locations),'ro');
    plot(x(1:index_DGPS_shuts_off_first,1),x_increment(1:index_DGPS_shuts_off_first,1),'bd','MarkerSize',4,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0.6 0.6 1]); 
    plot(x(index_DGPS_turns_on_last:end,1),x_increment(index_DGPS_turns_on_last:end,1),'bd','MarkerSize',4,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0.6 0.6 1]);
    plot(x,x_increment_pred,'m');
    plot(x_clean,x_increment_pred,'g');
    xlabel('x (m)');
    ylabel('x increments (m)');
    title('Skew removed from template.');
    legend('Raw data','Locations of jumps','DGPS lock at start','DGPS lock at end','Predicted data','Predicted with bias removed');
    
    
    figure(8824); clf;  hold on; grid minor;
    plot(x_clean,x_increments_clean,'r');
    plot(x_clean,x_increment_pred,'b');
    xlabel('x (m)');
    ylabel('x increments (m)');
    title('Skew and outliers removed from increments (clean).');
    legend('Increments with outliers and skew removed','Predicted Increments with skew removed');
    
    fcn_plotAxesLinkedTogetherByField;
end
%% Tell user we are leaving
if flag_do_debug
    % Show what we are doing
    fprintf(1,'\t\tExiting subfunction: %s\n',namestr);    
end

return

function [index_jump_locations,jump_mag] = fcn_calculateJumpLocations(x,x_increment_predicted)
flag_make_plots = 0;
jump_limit = 0.1;
N_standard_deviations = 3;

% calculate the biases and change in biases
x_increments = [0; diff(x(:,1))];
x_increments(1,1) = x_increment_predicted(1,1); % Replace the zero term
biases = cumsum(x_increments) - cumsum(x_increment_predicted);
change_in_biases = [0;diff(biases)];

% The standard deviation will be grossly skewed by jumps, so calculate a
% "clean" standard deviationt that does not have the jumps within the
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
    plot(time,jump_mag,'g');
    xlabel('Time (sec)');
    ylabel('Bias values (m)');
    legend('Biases','Locations of Jumps in Bias','Biases with jumps estimated and removed','Bias due to jumps only');
    title('Plot of biases');
end
return

