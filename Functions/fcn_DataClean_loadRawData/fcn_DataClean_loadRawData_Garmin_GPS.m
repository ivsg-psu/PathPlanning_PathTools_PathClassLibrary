function GPS_Garmin = fcn_DataClean_loadRawData_Garmin_GPS(d,data_source,flag_do_debug)

% This function is used to load the raw data collected with the Penn State Mapping Van.
% This is the GPS_Novatel data
% Input Variables:
%      d = raw data from GPS_Garmin(format:struct)
%      data_source = the data source of the raw data, can be 'mat_file' or 'database'(format:struct)
% Returned Results:
%      GPS_Garmin
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
    GPS_Garmin.ROS_Time       = d.Time';
    GPS_Garmin.centiSeconds  = 4; % This is sampled every 4 ms

    GPS_Garmin.Npoints        = length(GPS_Garmin.ROS_Time(:,1));
    deltaTSample              = mean(diff(GPS_Garmin.ROS_Time));
    GPS_Garmin.EmptyVector    = fcn_DataClean_fillEmptyStructureVector(GPS_Garmin); % Fill in empty vector (this is useful later)

    %GPS_Garmin.GPS_Time       = d.Time'*NaN;
    GPS_Garmin.Latitude       = d.Latitude';
    GPS_Garmin.Longitude      = d.Longitude';
    GPS_Garmin.Altitude       = d.Height';
    GPS_Garmin.xEast          = d.xEast';
    GPS_Garmin.yNorth         = d.yNorth';
    GPS_Garmin.zUp            = d.zUp';
    GPS_Garmin.velNorth       = [0; diff(GPS_Garmin.yNorth)]/deltaTSample;
    GPS_Garmin.velEast        = [0; diff(GPS_Garmin.xEast)]/deltaTSample;
    GPS_Garmin.velUp          = [0; diff(GPS_Garmin.zUp)]/deltaTSample;
    GPS_Garmin.velMagnitude   = sqrt(GPS_Garmin.velNorth.^2+GPS_Garmin.velEast.^2);
    % GPS_Garmin.numSatellites  = GPS_Garmin.EmptyVector;
    % GPS_Garmin.navMode        = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Roll_deg       = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Pitch_deg      = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Yaw_deg        = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Yaw_deg_Sigma  = GPS_Garmin.EmptyVector;

elseif strcmp(data_source,'database')
    GPS_Garmin.ROS_Time       = d.time;
    GPS_Garmin.centiSeconds  = 4; % This is sampled every 4 ms
    
    GPS_Garmin.Npoints        = length(GPS_Garmin.ROS_Time(:,1));
    deltaTSample              = mean(diff(GPS_Garmin.ROS_Time));
    GPS_Garmin.EmptyVector    = fcn_DataClean_fillEmptyStructureVector(GPS_Garmin); % Fill in empty vector (this is useful later)
    
    %GPS_Garmin.GPS_Time       = d.Time'*NaN;
    GPS_Garmin.Latitude       = d.latitude;
    GPS_Garmin.Longitude      = d.longitude;
    GPS_Garmin.Altitude       = d.altitude;
    GPS_Garmin.xEast          = d.xEast;
    GPS_Garmin.yNorth         = d.yNorth;
    GPS_Garmin.zUp            = d.zUp;
    GPS_Garmin.velNorth       = [0; diff(GPS_Garmin.yNorth)]/deltaTSample;
    GPS_Garmin.velEast        = [0; diff(GPS_Garmin.xEast)]/deltaTSample;
    GPS_Garmin.velUp          = [0; diff(GPS_Garmin.zUp)]/deltaTSample;
    GPS_Garmin.velMagnitude   = sqrt(GPS_Garmin.velNorth.^2+GPS_Garmin.velEast.^2);
    % GPS_Garmin.numSatellites  = GPS_Garmin.EmptyVector;
    % GPS_Garmin.navMode        = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Roll_deg       = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Pitch_deg      = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Yaw_deg        = GPS_Garmin.EmptyVector;
    % GPS_Garmin.Yaw_deg_Sigma  = GPS_Garmin.EmptyVector;
else
    error('Please indicate the data source')
end

GPS_Garmin.OneSigmaPos    = 3; % Typical 1-sigma in position, in meters

GPS_Garmin.xEast_Sigma    = GPS_Garmin.OneSigmaPos;
GPS_Garmin.yNorth_Sigma   = GPS_Garmin.OneSigmaPos;
GPS_Garmin.zUp_Sigma      = GPS_Garmin.OneSigmaPos;

% Calculate the increments
GPS_Garmin = fcn_DataClean_fillPositionIncrementsFromGPSPosition(GPS_Garmin);

% Calculate the apparent yaw angle from ENU positions
GPS_Garmin = fcn_DataClean_fillRawYawEstimatesFromGPSPosition(GPS_Garmin);

% Estimate the variance associated with the estimated yaw (Note: the Garmin
% is about 300 times less accurate than the Hemisphere, hence the 300
% factor.
GPS_Garmin.Yaw_deg_from_position_Sigma = fcn_DataClean_predictYawSigmaFromPosition(GPS_Garmin.xy_increments/300);

% Calculate the apparent yaw angle from ENU velocities
GPS_Garmin = fcn_DataClean_fillRawYawEstimatesFromGPSVelocity(GPS_Garmin);

% Estimate the variance associated with the estimated yaw
GPS_Garmin.Yaw_deg_from_velocity_Sigma = 100*ones(GPS_Garmin.Npoints,1);

% Close out the loading process
if flag_do_debug
    % Show what we are doing
    % Grab function name
    st = dbstack;
    namestr = st.name;
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return
