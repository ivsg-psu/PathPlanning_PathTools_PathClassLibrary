function cleaned = fcn_cleanIMUData(d)

%% First, copy all the data over
cleaned = d; 

%% Calculate time vector dat and transfer the standard items
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
    % Calculate the deltaT that the data should have, in GPS time (from
    % trigger)
    cleaned.GPS_Time_deltaT_target = 0.01*round(100*mean(diff(d.GPS_Time)));

end

% If both times exist, calculate difference
if isfield(d,'GPS_Time') && isfield(d,'ROS_Time')
    time_offsets = cleaned.ROS_Time - cleaned.GPS_Time;
    cleaned.ROS_to_GPS_Time_Offsets = time_offsets;
    cleaned.ROS_to_GPS_Time_Offsets_Sigma = std(time_offsets);
end

% %% Clean up the gyro measurements
% cleaned.ZGyro                  = d.ZGyro;
% cleaned.ZGyro_Sigma            = d.ZGyro_Sigma;
% 
% % STEP 1: Remove outliers
% cleaned.ZGyro = fcn_medianFilterData(cleaned.ZGyro,2*cleaned.ZGyro_Sigma);
% 
% % STEP 4: calculate the variance in the data
% cleaned.ZGyro_Sigma = std(diff(cleaned.ZGyro));
% 
% %% Clean up the XAccel measurements
% cleaned.XAccel                 = d.XAccel;
% cleaned.XAccel_Sigma           = d.XAccel_Sigma;
% cleaned.XAccel = fcn_medianFilterData(cleaned.XAccel,2*cleaned.XAccel_Sigma);
% cleaned.XAccel_Sigma = std(diff(cleaned.XAccel)); % Calculate the variance in the data
% 
% %% Clean up the YAccel measurements
% cleaned.YAccel                 = d.YAccel;
% cleaned.YAccel_Sigma           = d.YAccel_Sigma;
% cleaned.YAccel = fcn_medianFilterData(cleaned.YAccel,2*cleaned.YAccel_Sigma);
% cleaned.YAccel_Sigma = std(diff(cleaned.YAccel)); % Calculate the variance in the data
% 
% %% Clean up the ZAccel measurements
% cleaned.ZAccel                 = d.ZAccel;
% cleaned.ZAccel_Sigma           = d.ZAccel_Sigma;
% cleaned.ZAccel = fcn_medianFilterData(cleaned.ZAccel,2*cleaned.ZAccel_Sigma);
% cleaned.ZAccel_Sigma = std(diff(cleaned.ZAccel)); % Calculate the variance in the data

end

%% Subfunctions start here
function cleaned = fcn_medianFilterData(data,sigmas)
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
end

function unwrapped_angle = fcn_unwrapAngles(wrapped)
initial_angle = wrapped(1,1);
change_in_angle = [0; diff(wrapped)];
index_jumps = find(change_in_angle>180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)-360;
index_jumps = find(change_in_angle<-180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)+360;
unwrapped_angle = cumsum(change_in_angle) + initial_angle;

% Shift all data up or down 
mean_angle = mean(unwrapped_angle);
good_mean = mod(mean_angle,360);
shift = mean_angle - good_mean;
unwrapped_angle = unwrapped_angle - shift;
end

%% This function makes sure there are no missing samples...
function d = fcn_uniformSample(d)

initial_angle = wrapped(1,1);
change_in_angle = [0; diff(wrapped)];
index_jumps = find(change_in_angle>180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)-360;
index_jumps = find(change_in_angle<-180);
change_in_angle(index_jumps) = change_in_angle(index_jumps)+360;
unwrapped_angle = cumsum(change_in_angle) + initial_angle;

% Shift all data up or down 
mean_angle = mean(unwrapped_angle);
good_mean = mod(mean_angle,360);
shift = mean_angle - good_mean;
unwrapped_angle = unwrapped_angle - shift;
end
