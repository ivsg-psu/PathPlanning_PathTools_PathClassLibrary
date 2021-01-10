function Hemisphere = fcn_DataClean_loadRawData_Hemisphere(d,data_source,flag_do_debug)

% This function is used to load the raw data collected with the Penn State Mapping Van.
% This is the Hemisphere data
% Input Variables:
%      d = raw data from Hemisphere DGPS(format:struct)
%      data_source = the data source of the raw data, can be 'mat_file' or 'database'(format:struct)
% Returned Results:
%      Hemisphere
% Author: Liming Gao
% Created Date: 2020_11_15
% Modify Date: 2019_11_22
%
% Updates:
%
% To do lists:
% 1. check if it is reasonable for the calcualtion of Hemisphere.velMagnitude_Sigma
% 
%%
flag_plot = 0;

% the field name from mat_file is different from database, so we process
% them seperately
if strcmp(data_source,'mat_file')
    Hemisphere.ROS_Time         = d.Time';
    Hemisphere.GPS_Time         = d.GPSTimeOfWeek';
    Hemisphere.centiSeconds     = 5; % This is sampled every 5 ms
    
    Hemisphere.Npoints          = length(Hemisphere.ROS_Time(:,1));
    Hemisphere.EmptyVector      = fcn_DataClean_fillEmptyStructureVector(Hemisphere); % Fill in empty vector (this is useful later)
    
    Hemisphere.Latitude         = d.Latitude';
    Hemisphere.Longitude        = d.Longitude';
    Hemisphere.Altitude         = d.Height';
    Hemisphere.xEast            = d.xEast';
    Hemisphere.yNorth           = d.yNorth';
    Hemisphere.zUp              = d.zUp';
    Hemisphere.velNorth         = d.VNorth';
    Hemisphere.velEast          = d.VEast';
    Hemisphere.velUp            = d.VUp';
    Hemisphere.velMagnitude     = sqrt(Hemisphere.velNorth.^2 + Hemisphere.velEast.^2 + Hemisphere.velUp.^2);
    % for debugging - shows that the Hemisphere's velocity signal is horribly bad
    % figure;plot(Hemisphere.ROS_Time-Hemisphere.ROS_Time(1,1),Hemisphere.velMagnitude);
    % figure;plot(Hemisphere.ROS_Time-Hemisphere.ROS_Time(1,1),Hemisphere.velNorth);
    % figure;plot(Hemisphere.ROS_Time-Hemisphere.ROS_Time(1,1),Hemisphere.velEast);
    % figure;plot(Hemisphere.ROS_Time-Hemisphere.ROS_Time(1,1),Hemisphere.velUp);
    
    Hemisphere.velMagnitude_Sigma = std(Hemisphere.velMagnitude)*ones(length(Hemisphere.velMagnitude(:,1)),1);
    %Hemisphere.numSatellites    = Hemisphere.EmptyVector;
    Hemisphere.DGPS_is_active   = 1.00*(d.NavMode==6)';
    %Hemisphere.Roll_deg         = Hemisphere.EmptyVector;
    %Hemisphere.Pitch_deg        = Hemisphere.EmptyVector;
    %Hemisphere.Yaw_deg          = Hemisphere.EmptyVector;
    %Hemisphere.Yaw_deg_Sigma    = Hemisphere.EmptyVector;
    Hemisphere.OneSigmaPos      = d.StdDevResid';
    
elseif strcmp(data_source,'database')
    Hemisphere.ROS_Time         = d.time;
    Hemisphere.GPS_Time         = d.gps_seconds;
    Hemisphere.centiSeconds     = 5; % This is sampled every 5 ms
    
    Hemisphere.Npoints          = length(Hemisphere.ROS_Time(:,1));
    Hemisphere.EmptyVector      = fcn_DataClean_fillEmptyStructureVector(Hemisphere); % Fill in empty vector (this is useful later)
    
    Hemisphere.Latitude         = d.latitude;
    Hemisphere.Longitude        = d.longitude;
    Hemisphere.Altitude         = d.altitude;
    Hemisphere.xEast            = d.xeast;
    Hemisphere.yNorth           = d.ynorth;
    Hemisphere.zUp              = d.zup;
    Hemisphere.velNorth         = d.vnorth;
    Hemisphere.velEast          = d.veast;
    Hemisphere.velUp            = d.vup;
    Hemisphere.velMagnitude     = sqrt(Hemisphere.velNorth.^2 + Hemisphere.velEast.^2 + Hemisphere.velUp.^2);
    
    Hemisphere.velMagnitude_Sigma = std(Hemisphere.velMagnitude)*ones(length(Hemisphere.velMagnitude(:,1)),1);
    %Hemisphere.numSatellites    = Hemisphere.EmptyVector;
    Hemisphere.DGPS_is_active   = 1.00*(d.navmode==6);
    %Hemisphere.Roll_deg         = Hemisphere.EmptyVector;
    %Hemisphere.Pitch_deg        = Hemisphere.EmptyVector;
    %Hemisphere.Yaw_deg          = Hemisphere.EmptyVector;
    %Hemisphere.Yaw_deg_Sigma    = Hemisphere.EmptyVector;
    Hemisphere.OneSigmaPos      = d.stddevresid;
else
    error('Please indicate the data source')
end

% Update the variances in the position information, based on GPS status?
if 1==0
    Index_DGPS_active = Hemisphere.DGPS_is_active==1;
    Hemisphere.OneSigmaPos = 0.5 * ones(length(Hemisphere.EmptyVector),1); % Sigma value with DGPS inactive is roughly 50 cm, assuming WAAS
    Hemisphere.OneSigmaPos(Index_DGPS_active)    = 0.015; % Sigma value with DGPS active is roughly 1.5 cm.
end

Hemisphere.xEast_Sigma    = Hemisphere.OneSigmaPos;
Hemisphere.yNorth_Sigma   = Hemisphere.OneSigmaPos;
Hemisphere.zUp_Sigma      = Hemisphere.OneSigmaPos;

% Calculate the increments
Hemisphere = fcn_DataClean_fillPositionIncrementsFromGPSPosition(Hemisphere);

% Calculate the apparent yaw angle from ENU positions
Hemisphere = fcn_DataClean_fillRawYawEstimatesFromGPSPosition(Hemisphere);

% Estimate the variance associated with the estimated yaw
Hemisphere.Yaw_deg_from_position_Sigma = fcn_DataClean_predictYawSigmaFromPosition(Hemisphere.xy_increments);

% Estimate the variance associated with the estimated yaw based on velocity
Hemisphere.Yaw_deg_from_position_Sigma = 0*Hemisphere.xy_increments;  % Initialize array to zero
good_indices = find(Hemisphere.DGPS_is_active==1);
Hemisphere.Yaw_deg_from_position_Sigma(good_indices) = ...
    fcn_DataClean_predictYawSigmaFromPosition(Hemisphere.xy_increments(good_indices))/5;
bad_indices = find(Hemisphere.DGPS_is_active==0);
Hemisphere.Yaw_deg_from_position_Sigma(bad_indices)...
    = fcn_DataClean_predictYawSigmaFromPosition(Hemisphere.xy_increments(bad_indices));

% Calculate the apparent yaw angle from ENU velocities
Hemisphere = fcn_DataClean_fillRawYawEstimatesFromGPSVelocity(Hemisphere);


% Estimate the variance associated with the estimated yaw based on velocity
speeds = (Hemisphere.velNorth.^2+Hemisphere.velEast.^2).^0.5;
filt_speeds = medfilt1(speeds,20,'truncate');
Hemisphere.Yaw_deg_from_velocity_Sigma...
    = fcn_DataClean_predictYawSigmaFromVelocity(filt_speeds,Hemisphere.OneSigmaPos);

% The Hemisphere's yaw covariance grows much faster than position
% covariance, so we have to update this here to allow that...
% updates. See:
%     figure(454); plot(Hemisphere.GPS_Time - ceil(Hemisphere.GPS_Time(1,1)), Hemisphere.OneSigmaPos); xlim([539 546]);
% figure(46464);
% clf;
% plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1), Hemisphere.Yaw_deg_from_velocity,'r');
% hold on;
% plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1), Hemisphere.Yaw_deg_from_velocity+2*Hemisphere.Yaw_deg_from_velocity_Sigma,'Color',[1 0.5 0.5]);
% plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1), Hemisphere.Yaw_deg_from_velocity-2*Hemisphere.Yaw_deg_from_velocity_Sigma,'Color',[1 0.5 0.5]);
% xlim([539 546]);

% This should show the variance inflate inbetween updates. Because it is
% not doing this correctly, we'll need to estimate this and fix it..
time_offset = -0.35;
time = Hemisphere.GPS_Time  + time_offset;
t_fraction = time - floor(time);
%figure(46346); plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1),t_fraction); xlim([539 546]);

a = 0.03;  % Units are centimeters - this is an estimate of how much the GPS will linearly drift after 1 second, without corrections
b = .7;  % Units are centimeters per t^0.5. This reprents the random walk rate of the GPS without corrections
new_variance = Hemisphere.OneSigmaPos + a*t_fraction + b*t_fraction.^0.5;

% Recalculate based on quadratic growth model given above
Hemisphere.Yaw_deg_from_velocity_Sigma...
    = fcn_DataClean_predictYawSigmaFromVelocity(filt_speeds,new_variance);
if flag_plot
    figure(46464);
    clf;
    plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1), Hemisphere.Yaw_deg_from_velocity,'r');
    hold on;
    plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1), Hemisphere.Yaw_deg_from_velocity+2*Hemisphere.Yaw_deg_from_velocity_Sigma,'Color',[1 0.5 0.5]);
    plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1), Hemisphere.Yaw_deg_from_velocity-2*Hemisphere.Yaw_deg_from_velocity_Sigma,'Color',[1 0.5 0.5]);
    xlim([539 546]);
end

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