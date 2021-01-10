function GPS_Novatel = fcn_DataClean_loadRawData_Novatel_GPS(d,Hemisphere,data_source,flag_do_debug)

% This function is used to load the raw data collected with the Penn State Mapping Van.
% This is the GPS_Novatel data
% Input Variables:
%      d = raw data from GPS_Novatel(format:struct)
%      Hemisphere = the data from Hemisphere GPS, used to estimate the
%                   GPS_Novatel sigma (format:struct)
%      data_source = the data source of the raw data, can be 'mat_file' or 'database'(format:struct)
% Returned Results:
%      GPS_Novatel
% Author: Liming Gao
% Created Date: 2020_12_07
%
%
% Updates:
%
% To do lists:
% 1. check if it is reasonable to select data from the second d.Time(2:end)';
% 2. check the Yaw_deg between matfile and database
% 3. Hemisphere = d_out;  %%update the interpolated values to raw data?
%%
% the field name from mat_file is different from database, so we process
% them seperately
if strcmp(data_source,'mat_file')
    % % Note: the Novatel and Hemisphere are almost perfectly time aligned, if
    % % dropping the first data point in Novatel (uncomment the following to
    % % see)
    % Hemisphere.GPS_Time(1,1)
    % % ans =
    % %           242007.249999977
    % d.Seconds(1,2)
    % % ans =
    % %              242007.248687
    % % This is why all the vectors below start at 2, not 1
    GPS_Novatel.ROS_Time       = d.Time(2:end)';
    GPS_Novatel.GPS_Time       = d.Seconds(2:end)';
    GPS_Novatel.centiSeconds   = 5; % This is sampled every 5 ms
    
    GPS_Novatel.Npoints        = length(GPS_Novatel.ROS_Time(:,1));
    GPS_Novatel.EmptyVector    = fcn_DataClean_fillEmptyStructureVector(GPS_Novatel); % Fill in empty vector (this is useful later)
    
    GPS_Novatel.Latitude       = d.Latitude(2:end)';
    GPS_Novatel.Longitude      = d.Longitude(2:end)';
    GPS_Novatel.Altitude       = d.Height(2:end)';
    GPS_Novatel.xEast          = d.xEast(2:end)';
    GPS_Novatel.yNorth         = d.yNorth(2:end)';
    GPS_Novatel.zUp            = d.zUp(2:end)';
    GPS_Novatel.velNorth       = d.NorthVelocity(2:end)';
    GPS_Novatel.velEast        = d.EastVelocity(2:end)';
    GPS_Novatel.velUp          = d.UpVelocity(2:end)';
    GPS_Novatel.velMagnitude   = sqrt(d.NorthVelocity(2:end)'.^2+d.EastVelocity(2:end)'.^2);
    GPS_Novatel.velMagnitude_Sigma = std(diff(GPS_Novatel.velMagnitude))*ones(length(GPS_Novatel.velMagnitude(:,1)),1);
    GPS_Novatel.DGPS_is_active = zeros(GPS_Novatel.Npoints,1);
    GPS_Novatel.numSatellites  = GPS_Novatel.EmptyVector;
    GPS_Novatel.navMode        = GPS_Novatel.EmptyVector;
    GPS_Novatel.Roll_deg       = d.Roll(2:end)';
    GPS_Novatel.Pitch_deg      = d.Pitch(2:end)';
    GPS_Novatel.Yaw_deg        = -d.Azimuth(2:end)'+360+90; % Notice sign flip and phase shift due to coord convention and mounting

elseif strcmp(data_source,'database')
    GPS_Novatel.ROS_Time       = d.time;
    GPS_Novatel.GPS_Time       = d.gps_seconds;
    GPS_Novatel.centiSeconds   = 5; % This is sampled every 5 ms
    
    GPS_Novatel.Npoints        = length(GPS_Novatel.ROS_Time(:,1));
    GPS_Novatel.EmptyVector    = fcn_DataClean_fillEmptyStructureVector(GPS_Novatel); % Fill in empty vector (this is useful later)
    
    GPS_Novatel.Latitude       = d.latitude;
    GPS_Novatel.Longitude      = d.longitude;
    GPS_Novatel.Altitude       = d.altitude;
    GPS_Novatel.xEast          = d.xEast;
    GPS_Novatel.yNorth         = d.yNorth;
    GPS_Novatel.zUp            = d.zUp;
    GPS_Novatel.velNorth       = d.north_velocity;
    GPS_Novatel.velEast        = d.east_velocity;
    GPS_Novatel.velUp          = d.up_velocity;
    GPS_Novatel.velMagnitude   = sqrt(d.north_velocity.^2 + d.east_velocity.^2);
    GPS_Novatel.velMagnitude_Sigma = std(diff(GPS_Novatel.velMagnitude))*ones(length(GPS_Novatel.velMagnitude(:,1)),1);
    GPS_Novatel.DGPS_is_active = zeros(GPS_Novatel.Npoints,1);
    GPS_Novatel.numSatellites  = GPS_Novatel.EmptyVector;
    GPS_Novatel.navMode        = GPS_Novatel.EmptyVector;
    GPS_Novatel.Roll_deg       = rad2deg(d.roll);
    GPS_Novatel.Pitch_deg      = rad2deg(d.pitch);
    GPS_Novatel.Yaw_deg        = rad2deg(d.yaw)+360; % in database the yaw is 90 - yaw
else
    error('Please indicate the data source')
end

GPS_Novatel.Yaw_deg_Sigma  = 0.1 * ones(GPS_Novatel.Npoints,1); % Units are deg. Constant due to IMU on Novatel
GPS_Novatel.OneSigmaPos    = 0.30; % Typical 1-sigma in position, in meters

% The sigma values for the ENU coordinates are not given. But the
% Hemisphere is known nearly precisely due to differential corrections
% (usually). So we define the Novatel's variance as the error between it
% and the Hemisphere, plus a multiple times the Hemisphere's variance.

% Before calculating the sigma, interpolate  hemisphere data using
% GPS_Novatel time sequence to ensure they have the same data length
d_in = Hemisphere; % Create a temporary data structure
names = fieldnames(d_in);

for i_data = 1: length(names)
    
    data_name = names{i_data};
    d_mid = d_in.(data_name); % extract field  data
    if length(d_mid) == length(d_in.ROS_Time)
        d_out.(data_name) = interp1(d_in.ROS_Time, d_mid, GPS_Novatel.GPS_Time, 'linear', 'extrap');
    else
        d_out.(data_name) = d_mid;
    end
end
Hemisphere = d_out;  %%update this to raw data?

multiple = 3;
GPS_Novatel.xEast_Sigma    = abs(GPS_Novatel.xEast - Hemisphere.xEast) + multiple * Hemisphere.xEast_Sigma;
GPS_Novatel.yNorth_Sigma   = abs(GPS_Novatel.yNorth - Hemisphere.yNorth) + multiple * Hemisphere.yNorth_Sigma;
GPS_Novatel.zUp_Sigma      = abs(GPS_Novatel.zUp - Hemisphere.zUp) + multiple * Hemisphere.zUp_Sigma;

% Below is old method:
% scaling = 0.01;
% GPS_Novatel.xEast_Sigma    = GPS_Novatel.OneSigmaPos * scaling;
% GPS_Novatel.yNorth_Sigma   = GPS_Novatel.OneSigmaPos * scaling;
% GPS_Novatel.zUp_Sigma      = GPS_Novatel.OneSigmaPos * scaling;

% Calculate the increments
GPS_Novatel = fcn_DataClean_fillPositionIncrementsFromGPSPosition(GPS_Novatel);

% Calculate the apparent yaw angle from ENU positions
GPS_Novatel = fcn_DataClean_fillRawYawEstimatesFromGPSPosition(GPS_Novatel);

% Estimate the variance associated with the estimated yaw
GPS_Novatel.Yaw_deg_from_position_Sigma = 0.5*ones(GPS_Novatel.Npoints,1);


% Calculate the apparent yaw angle from ENU velocities
GPS_Novatel = fcn_DataClean_fillRawYawEstimatesFromGPSVelocity(GPS_Novatel);

% Estimate the variance associated with the estimated yaw
speeds = (GPS_Novatel.velNorth.^2+GPS_Novatel.velEast.^2).^0.5;
filt_speeds = medfilt1(speeds,20,'truncate');
GPS_Novatel.Yaw_deg_from_velocity_Sigma...
    = fcn_DataClean_predictYawSigmaFromVelocity(filt_speeds,GPS_Novatel.OneSigmaPos);

clear d %clear temp variable


% Close out the loading process
if flag_do_debug
    % Show what we are doing
    % Grab function name
    st = dbstack;
    namestr = st.name;
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return
