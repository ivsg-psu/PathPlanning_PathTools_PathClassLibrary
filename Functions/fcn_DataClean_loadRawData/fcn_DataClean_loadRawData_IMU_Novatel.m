function IMU_Novatel = fcn_DataClean_loadRawData_IMU_Novatel(d,data_source,flag_do_debug)

% This function is used to load the raw data collected with the Penn State Mapping Van.
% This is the GPS_Novatel data
% Input Variables:
%      d = raw data from GPS_Novatel(format:struct)
%      data_source = the data source of the raw data, can be 'mat_file' or 'database'(format:struct)
%
% Returned Results:
%      IMU_Novatel
% Author: Liming Gao
% Created Date: 2020_12_07
%
%
% Updates:
%
% To do lists:
% 1.
%
%

%%
% the field name from mat_file is different from database, so we process
% them seperately
if strcmp(data_source,'mat_file')
    IMU_Novatel.ROS_Time      = d.Time';
    IMU_Novatel.centiSeconds  = 1; % This is sampled every 1 ms
    
    IMU_Novatel.Npoints       = length(IMU_Novatel.ROS_Time(:,1));
    IMU_Novatel.EmptyVector   = fcn_DataClean_fillEmptyStructureVector(IMU_Novatel); % Fill in empty vector (this is useful later)
    IMU_Novatel.GPS_Time      = d.Seconds';
    IMU_Novatel.deltaT_ROS    = mean(diff(IMU_Novatel.ROS_Time));
    IMU_Novatel.deltaT_GPS    = mean(diff(IMU_Novatel.GPS_Time));
    IMU_Novatel.IMUStatus     = d.IMUStatus';
    IMU_Novatel.XAccel        = d.XAccel';
    IMU_Novatel.YAccel        = d.YAccel';
    IMU_Novatel.ZAccel        = d.ZAccel';
    IMU_Novatel.XGyro         = d.XGyro';
    IMU_Novatel.YGyro         = d.YGyro';
    IMU_Novatel.ZGyro         = d.ZGyro';
elseif strcmp(data_source,'database')
    IMU_Novatel.ROS_Time      = d.time;
    IMU_Novatel.centiSeconds  = 1; % This is sampled every 1 ms
    
    IMU_Novatel.Npoints       = length(IMU_Novatel.ROS_Time(:,1));
    IMU_Novatel.EmptyVector   = fcn_DataClean_fillEmptyStructureVector(IMU_Novatel); % Fill in empty vector (this is useful later)
    IMU_Novatel.GPS_Time      = d.gps_seconds;
    IMU_Novatel.deltaT_ROS    = mean(diff(IMU_Novatel.ROS_Time));
    IMU_Novatel.deltaT_GPS    = mean(diff(IMU_Novatel.GPS_Time));
    IMU_Novatel.IMUStatus     = d.status;
    IMU_Novatel.XAccel        = d.x_acceleration;
    IMU_Novatel.YAccel        = d.y_acceleration;
    IMU_Novatel.ZAccel        = d.z_acceleration;
    IMU_Novatel.XGyro         = d.x_angular_velocity;
    IMU_Novatel.YGyro         = d.y_angular_velocity;
    IMU_Novatel.ZGyro         = d.z_angular_velocity;
else
    error('Please indicate the data source')
end

% Close out the loading process
if flag_do_debug
    % Show what we are doing
    % Grab function name
    st = dbstack;
    namestr = st.name;
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return
