function cleaned = fcn_DataClean_cleanGPSData(d)

flag_do_debug = 0;

%% First, copy all the data over
cleaned = d;

%% TIME CALCULATIONS: Calculate time vector data and transfer the standard items
if isfield(d,'ROS_Time')
    cleaned.ROS_Time               = d.ROS_Time;
    cleaned.ROS_Time_Sigma         = std(diff(d.ROS_Time));
    cleaned.ROS_Time_deltaT        = mean(diff(d.ROS_Time));
    cleaned.ROS_Time_deltaT_target = 0.01*round(100*cleaned.ROS_Time_deltaT);  % Calculate the deltaT that the data should have (from trigger)
end

if isfield(d,'GPS_Time')
    cleaned.GPS_Time               = d.GPS_Time;
    cleaned.GPS_Time_Sigma         = std(diff(d.GPS_Time));
    cleaned.GPS_Time_deltaT        = mean(diff(d.GPS_Time));
    cleaned.GPS_Time_deltaT_target = 0.01*round(100*cleaned.GPS_Time_deltaT);  % Calculate the deltaT that the data should have (from trigger)
end

cleaned.EmptyVector            = d.EmptyVector;



%% YAW POSITION CALCULATIONS: Clean the position-based Yaw values

t = d.ROS_Time;

%% STEP 1: Remove jumps caused by wrapping
raw_yaw = d.Yaw_deg_from_position;
unwrapped_yaw = fcn_DataClean_unwrapAngles(raw_yaw);

if 1==flag_do_debug
    figure(46473);
    clf;
    plot(t - t(1),raw_yaw,'r');
    hold on;
    plot(t - t(1),unwrapped_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Raw angle','Unwrapped');
    h_title = title('GPS_Hemisphere Yaw_deg_from_position');
    set(h_title, 'Interpreter', 'none');
end

%% STEP 2: Remove solid zeros
good_indices = find(abs(unwrapped_yaw)>0.000000001); % Set a very low tolerance...
dezeroed_yaw = fcn_DataClean_replaceBadIndicesWithNearestGood(unwrapped_yaw,good_indices);

if 1==flag_do_debug
    figure(9768675);
    clf;
    plot(t - t(1),unwrapped_yaw,'r');
    hold on;
    plot(t - t(1),dezeroed_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Unwrapped angle','Dezeroed angle');
    h_title = title('GPS_Hemisphere Yaw_deg_from_position');
    set(h_title, 'Interpreter', 'none');
end

%% STEP 3: Perform kinematic filtering
U = d.velMagnitude;  % This is the forward velocity
U_limit = 1; % This is 1 m/s, or about 2 mph
kinematic_yaw = fcn_DataClean_enforceKinematicYawLimits(dezeroed_yaw,t, U,U_limit);


if 1==flag_do_debug
    figure(565758);
    clf;
    plot(t - t(1),dezeroed_yaw,'r');
    hold on;
    plot(t - t(1),kinematic_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Dezeroed angle','Kinematic angle');
    h_title = title('GPS_Hemisphere Yaw_deg_from_position');
    set(h_title, 'Interpreter', 'none');
end

cleaned.Yaw_deg_from_position = kinematic_yaw;


%% STEP 4: calculate the variance in the data
cleaned.Yaw_deg_from_position_Sigma = d.Yaw_deg_from_position_Sigma;



%% Clean the velocity-based Yaw values

raw_yaw = d.Yaw_deg_from_velocity;


%% STEP 1: Remove situations where vehicle is stationary, 
% causing yaw problems
good_indices = find(d.velMagnitude>0.05); % There's enough motion to calculate yaw angle
nonstationary_yaw = fcn_DataClean_replaceBadIndicesWithNearestGood(raw_yaw,good_indices);


if 1==flag_do_debug
    figure(22565633);
    clf;
    plot(t - t(1),raw_yaw,'r');
    hold on;
    plot(t - t(1),nonstationary_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Raw angle','Nonstationary angle');
    h_title = title('GPS_Hemisphere Yaw_deg_from_velocity');
    set(h_title, 'Interpreter', 'none');
end

%% STEP 3: unwrap the yaw jumps out of the data
unwrapped_yaw = fcn_DataClean_unwrapAngles(nonstationary_yaw);


if 1==flag_do_debug
    figure(2246473);
    clf;
    plot(t - t(1),nonstationary_yaw,'r');
    hold on;
    plot(t - t(1),unwrapped_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Nonstationary','Unwrapped');
    h_title = title('GPS_Hemisphere Yaw_deg_from_velocity');
    set(h_title, 'Interpreter', 'none');
end


%% STEP 3: Perform kinematic filtering
U = d.velMagnitude;  % This is the forward velocity
U_limit = 0.05; % This is 0.05 m/s, or about 0.1 mph - note that it is much smaller here because velocity measurements are MUCH more accurate than position for yaw angle calculations
kinematic_yaw = fcn_DataClean_enforceKinematicYawLimits(unwrapped_yaw,t, U,U_limit);
kinematic_yaw = fcn_DataClean_unwrapAngles(kinematic_yaw);

if 1==flag_do_debug
    figure(22565758);
    clf;
    plot(t - t(1),unwrapped_yaw,'r');
    hold on;
    plot(t - t(1),kinematic_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Unwrapped angle','Kinematic angle');
    h_title = title('GPS_Hemisphere Yaw_deg_from_velocity');
    set(h_title, 'Interpreter', 'none');
end

cleaned.Yaw_deg_from_velocity = kinematic_yaw;

%% STEP 4: calculate the variance in the data
cleaned.Yaw_deg_from_velocity_Sigma = d.Yaw_deg_from_velocity_Sigma;

%% Fix jumps in position


%% STEP 2: Remove solid zeros
good_indices = find(abs(unwrapped_yaw)>0.000000001); % Set a very low tolerance...
dezeroed_yaw = fcn_DataClean_replaceBadIndicesWithNearestGood(unwrapped_yaw,good_indices);

if 1==flag_do_debug
    figure(9768675);
    clf;
    plot(t - t(1),unwrapped_yaw,'r');
    hold on;
    plot(t - t(1),dezeroed_yaw,'k');
    xlabel('ROS Time (sec)');
    ylabel('Yaw angle (deg)');
    legend('Unwrapped angle','Dezeroed angle');
    h_title = title('GPS_Hemisphere Yaw_deg_from_position');
    set(h_title, 'Interpreter', 'none');
end

return

%%
function cleaned = fcn_DataClean_medianFilterData(data,sigmas) %#ok<DEFNU>
data_median = medfilt1(data,7,'truncate');
sigma_median = medfilt1(sigmas,7,'truncate');

% Calculate bounds
highest_expected_data = data_median + sigma_median;
lowest_expected_data = data_median - sigma_median;

% For debugging:
% figure; plot(data,'b'); hold on; plot(data_median,'c'); plot(highest_expected_data,'r'); plot(lowest_expected_data,'r');

% Find outliers and remove them via the median filter
out_of_bounds = [...
    find(data>highest_expected_data);
    find(data<lowest_expected_data)];

cleaned = data;
cleaned(out_of_bounds) = data_median(out_of_bounds);

% For debugging:
% figure;  plot(cleaned,'b'); hold on; plot(highest_expected_data,'r'); plot(lowest_expected_data,'r');
return

%%
function unwrapped_angle = fcn_DataClean_unwrapAngles(wrapped)
% Sometimes this function is called with NaN values, which will give bad
% results, so we need to fill these in
do_debug = 0;
if 1== do_debug % For debugging
    figure(3737);
    clf;
    hold on;
    grid minor;
    plot(wrapped);
end
wrapped = fillmissing(wrapped,'previous');

if 1== do_debug % For debugging
    figure(3737);
    plot(wrapped,'r');
end



initial_angle = wrapped(1,1); % Grab the first angle
change_in_angle = [0; diff(wrapped)]; % Use diff to find jumps
index_jumps = find(change_in_angle>180); % Tag jumps greater than 180
change_in_angle(index_jumps) = change_in_angle(index_jumps)-360; % Shift these down
index_jumps = find(change_in_angle<-180); % Tag jumps less than -180
change_in_angle(index_jumps) = change_in_angle(index_jumps)+360; % Shift these up
unwrapped_angle = cumsum(change_in_angle) + initial_angle; % Re-add jumps

% After above process, data may be shifted up or down by multiples of 360,
% so shift all data back so mean is between zero and 360
mean_angle = mean(unwrapped_angle);
good_mean = mod(mean_angle,360);
shift = mean_angle - good_mean;
unwrapped_angle = unwrapped_angle - shift;
return

%%
function new_data = fcn_DataClean_replaceBadIndicesWithNearestGood(data,good_indices)

if ~isempty(good_indices)
    new_data = data; % Pre-fill new data
    all_indices = (1:length(data))'; % Grab all indices
    moved_enough = ismember(all_indices,good_indices); % Tag bad ones
    for i = 1:length(moved_enough) % Loop through indices
        if 0==moved_enough(i) % Check if bad one
            % Find closest index to stationary one that can work
            distances = (i-good_indices).^2;
            [~,index] = min(distances);
            best_index = good_indices(index);
            new_data(i,1) = data(best_index,1);
        end
    end
else
    % The good_indices vector is empty
    warning('Empty vector found for good indices - unclear how to fix.');
    new_data = data;
end
return

%% 
function kinematic_yaw = fcn_DataClean_enforceKinematicYawLimits(yaw,t, U,U_limit)
U(U<U_limit) = 0;  % Set the slow speeds to zero - can't trust GPS slower than this
L = 2.5; % Typical wheelbase in meters
maxsteer = 40*pi/180; % Maximum steering angle magnitude is about 40 degrees
maxyawrate = U * maxsteer/L; % The kinematics limit

kinematic_yaw = yaw; % Initialize the variable
for i_time = 2:length(kinematic_yaw) % Loop through variables
    delta_t = t(i_time) - t(i_time-1); % Calculate delta t
    yaw_max = kinematic_yaw(i_time-1,1) + maxyawrate(i_time,1)*delta_t*180/pi;
    yaw_min = kinematic_yaw(i_time-1,1) - maxyawrate(i_time,1)*delta_t*180/pi;
    
    if kinematic_yaw(i_time,1) > yaw_max
        kinematic_yaw(i_time,1) = yaw_max;
    elseif kinematic_yaw(i_time,1) < yaw_min
        kinematic_yaw(i_time,1) = yaw_min;
    end
end
return

