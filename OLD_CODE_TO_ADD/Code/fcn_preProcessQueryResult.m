function rawData = fcn_preProcessQueryResult(data_struct)

%{
flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Starting load of rawData structure from source files.\n');    
end

%% Load the data
data{1}.filename = filename;% 'Route_Wahba.mat';
data{1}.variable_names ={variable_names}; % {'Route_WahbaLoop'};

for i_data = 1:length(data)
    ith_filename = data{i_data}.filename;
    ith_variable_name = data{i_data}.variable_names{1}; % Need to do this as a loop if more than one variable
    if flag_do_debug
        % Show what we are doing
        fprintf(1,'Source file: %s is being used to load variable %s\n',ith_filename,ith_variable_name);
    end
    data_name = load(ith_filename,ith_variable_name);
end
    data_struct = data_name.( ith_variable_name); %Accessing Data Using Dynamic Field Names
%}
 %% Process data from the GPS_2019 mat file - This is the Hemisphere
d = data_struct.Hemisphere_DGPS; % Create a temporary data structure

Hemisphere.ROS_Time         = d.time;
Hemisphere.GPS_Time         = d.gps_seconds;
Hemisphere.centiSeconds     = 5; % This is sampled every 5 ms

Hemisphere.Npoints          = length(Hemisphere.ROS_Time(:,1));
Hemisphere.EmptyVector      = fcn_fillEmptyStructureVector(Hemisphere); % Fill in empty vector (this is useful later)

Hemisphere.Latitude         = d.latitude;
Hemisphere.Longitude        = d.longitude;
Hemisphere.Altitude         = d.altitude;
Hemisphere.xEast            = d.xEast;
Hemisphere.yNorth           = d.yNorth;
Hemisphere.zUp              = d.zUp;
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
Hemisphere.xEast_Sigma    = Hemisphere.OneSigmaPos;
Hemisphere.yNorth_Sigma   = Hemisphere.OneSigmaPos;
Hemisphere.zUp_Sigma      = Hemisphere.OneSigmaPos;

% Calculate the increments
Hemisphere = fcn_fillPositionIncrementsFromGPSPosition(Hemisphere);

% Calculate the apparent yaw angle from ENU positions
Hemisphere = fcn_fillRawYawEstimatesFromGPSPosition(Hemisphere);

% Estimate the variance associated with the estimated yaw
Hemisphere.Yaw_deg_from_position_Sigma = fcn_predictYawSigmaFromPosition(Hemisphere.xy_increments);

% Estimate the variance associated with the estimated yaw based on velocity
Hemisphere.Yaw_deg_from_position_Sigma = 0*Hemisphere.xy_increments;  % Initialize array to zero
good_indices = find(Hemisphere.DGPS_is_active==1);
Hemisphere.Yaw_deg_from_position_Sigma(good_indices) = fcn_predictYawSigmaFromPosition(Hemisphere.xy_increments(good_indices))/5;

bad_indices = find(Hemisphere.DGPS_is_active==0);
Hemisphere.Yaw_deg_from_position_Sigma(bad_indices) = fcn_predictYawSigmaFromPosition(Hemisphere.xy_increments(bad_indices));

% Calculate the apparent yaw angle from ENU velocities
Hemisphere = fcn_fillRawYawEstimatesFromGPSVelocity(Hemisphere);

% Estimate the variance associated with the estimated yaw based on velocity
speeds = (Hemisphere.velNorth.^2 + Hemisphere.velEast.^2).^0.5;
filt_speeds = medfilt1(speeds,20,'truncate');
Hemisphere.Yaw_deg_from_velocity_Sigma = fcn_predictYawSigmaFromVelocity(filt_speeds,Hemisphere.OneSigmaPos);

% This should show the variance inflate inbetween updates. Because it is
% not doing this correctly, we'll need to estimate this and fix it..
time_offset = -0.35;
time = Hemisphere.GPS_Time  + time_offset;
t_fraction = time - floor(time);

a = 0.03;  % Units are centimeters - this is an estimate of how much the GPS will linearly drift after 1 second, without corrections
b = .7;  % Units are centimeters per t^0.5. This reprents the random walk rate of the GPS without corrections
new_variance = Hemisphere.OneSigmaPos + a*t_fraction + b*t_fraction.^0.5;

% Recalculate based on quadratic growth model given above
Hemisphere.Yaw_deg_from_velocity_Sigma = fcn_predictYawSigmaFromVelocity(filt_speeds,new_variance);

rawData.GPS_Hemisphere = Hemisphere;

clear d %clear temp variable 
%% Process data from Novatel
            
d = data_struct.GPS_Novatel; % Create a temporary data structure

% Fill in the fields
GPS_Novatel.ROS_Time       = d.time;
GPS_Novatel.GPS_Time       = d.gps_seconds;
GPS_Novatel.centiSeconds   = 5; % This is sampled every 5 ms

GPS_Novatel.Npoints        = length(GPS_Novatel.ROS_Time(:,1));
GPS_Novatel.EmptyVector    = fcn_fillEmptyStructureVector(GPS_Novatel); % Fill in empty vector (this is useful later)

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
%GPS_Novatel.Yaw_deg        = -d.Azimuth(2:end)'+360+90; % Notice sign flip and phase shift due to coord convention and mounting
GPS_Novatel.Yaw_deg_Sigma  = 0.1 * ones(GPS_Novatel.Npoints,1); % Units are deg. Constant due to IMU on Novatel
GPS_Novatel.OneSigmaPos    = 0.30; % Typical 1-sigma in position, in meters

% The sigma values for the ENU coordinates are not given. But the
% Hemisphere is known nearly precisely due to differential corrections
% (usually). So we define the Novatel's variance as the error between it
% and the Hemisphere, plus a multiple times the Hemisphere's variance.

% Before calculating the sigma, interpolate GPS_Novatel data using
% hemisphere time sequence to ensure they have the same data length
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

% Calculate the increments
GPS_Novatel = fcn_fillPositionIncrementsFromGPSPosition(GPS_Novatel);

% Calculate the apparent yaw angle from ENU positions
GPS_Novatel = fcn_fillRawYawEstimatesFromGPSPosition(GPS_Novatel);

% Estimate the variance associated with the estimated yaw
GPS_Novatel.Yaw_deg_from_position_Sigma = 0.5*ones(GPS_Novatel.Npoints,1);

% Calculate the apparent yaw angle from ENU velocities
GPS_Novatel = fcn_fillRawYawEstimatesFromGPSVelocity(GPS_Novatel);

% Estimate the variance associated with the estimated yaw
speeds = (GPS_Novatel.velNorth.^2+GPS_Novatel.velEast.^2).^0.5;
filt_speeds = medfilt1(speeds,20,'truncate');
GPS_Novatel.Yaw_deg_from_velocity_Sigma = fcn_predictYawSigmaFromVelocity(filt_speeds,GPS_Novatel.OneSigmaPos);

rawData.GPS_Novatel = GPS_Novatel;
clear d_out d %clear temp variable 
%% Process data from the Garmin
d = data_struct.Garmin_GPS; % Create a temporary data structure

GPS_Garmin.ROS_Time       = d.time;
GPS_Garmin.centiSeconds  = 4; % This is sampled every 4 ms

GPS_Garmin.Npoints        = length(GPS_Garmin.ROS_Time(:,1));
deltaTSample              = mean(diff(GPS_Garmin.ROS_Time));
GPS_Garmin.EmptyVector    = fcn_fillEmptyStructureVector(GPS_Garmin); % Fill in empty vector (this is useful later)

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
GPS_Garmin.OneSigmaPos    = 3; % Typical 1-sigma in position, in meters

GPS_Garmin.xEast_Sigma    = GPS_Garmin.OneSigmaPos;
GPS_Garmin.yNorth_Sigma   = GPS_Garmin.OneSigmaPos;
GPS_Garmin.zUp_Sigma      = GPS_Garmin.OneSigmaPos;

% Calculate the increments
GPS_Garmin = fcn_fillPositionIncrementsFromGPSPosition(GPS_Garmin);

% Calculate the apparent yaw angle from ENU positions
GPS_Garmin = fcn_fillRawYawEstimatesFromGPSPosition(GPS_Garmin);

% Estimate the variance associated with the estimated yaw (Note: the Garmin
% is about 300 times less accurate than the Hemisphere, hence the 300
% factor.
GPS_Garmin.Yaw_deg_from_position_Sigma = fcn_predictYawSigmaFromPosition(GPS_Garmin.xy_increments/300);

% Calculate the apparent yaw angle from ENU velocities
GPS_Garmin = fcn_fillRawYawEstimatesFromGPSVelocity(GPS_Garmin);

% Estimate the variance associated with the estimated yaw
GPS_Garmin.Yaw_deg_from_velocity_Sigma = 100*ones(GPS_Garmin.Npoints,1);

rawData.GPS_Garmin = GPS_Garmin;

%% Process data from the Novatel IMU
d = data_struct.Novatel_IMU; % Create a temporary data structure
IMU_Novatel.ROS_Time      = d.time;
IMU_Novatel.centiSeconds  = 1; % This is sampled every 1 ms

IMU_Novatel.Npoints       = length(IMU_Novatel.ROS_Time(:,1));
IMU_Novatel.EmptyVector   = fcn_fillEmptyStructureVector(IMU_Novatel); % Fill in empty vector (this is useful later)
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

rawData.IMU_Novatel = IMU_Novatel;

%% Process data from the ADIS IMU
d = data_struct.adis_IMU_data; % Create a temporary data structure
IMU_ADIS.ROS_Time      = d.time;
IMU_ADIS.centiSeconds  = 1; % This is sampled every 1 ms
IMU_ADIS.Npoints       = length(IMU_ADIS.ROS_Time(:,1));
IMU_ADIS.EmptyVector   = fcn_fillEmptyStructureVector(IMU_ADIS); % Fill in empty vector (this is useful later)
%IMU_ADIS.GPS_Time      = IMU_ADIS.EmptyVector;
IMU_ADIS.deltaT_ROS    = mean(diff(IMU_ADIS.ROS_Time));
%IMU_ADIS.deltaT_GPS    = mean(diff(IMU_ADIS.GPS_Time));
%IMU_ADIS.IMUStatus     = IMU_ADIS.EmptyVector;

% x and y reversed? x -> -y?
IMU_ADIS.XAccel        = -d.y_acceleration;
IMU_ADIS.YAccel        = d.x_acceleration;
IMU_ADIS.ZAccel        = d.z_acceleration;

% x and y swapped: x -> y
IMU_ADIS.XGyro         = d.y_angular_velocity; 
IMU_ADIS.YGyro         = d.x_angular_velocity;
IMU_ADIS.ZGyro         = d.z_angular_velocity;

rawData.IMU_ADIS = IMU_ADIS;

%% Process data from the steering sensor - the sensor stinks, so we won't use it
d = data_struct.Steering_angle; % Create a temporary data structure

Input_Steering.ROS_Time        = d.time;
Input_Steering.centiSeconds    = 1; % This is sampled every 1 ms
Input_Steering.Npoints         = length(Input_Steering.ROS_Time(:,1));
Input_Steering.EmptyVector     = fcn_fillEmptyStructureVector(Input_Steering); % Fill in empty vector (this is useful later)
%Input_Steering.GPS_Time        = Input_Steering.EmptyVector;
Input_Steering.deltaT_ROS      = mean(diff(Input_Steering.ROS_Time));
%Input_Steering.deltaT_GPS      = mean(diff(Input_Steering.GPS_Time));
Input_Steering.LeftAngle       = -1*d.left_angle;
Input_Steering.RightAngle      = 1*d.right_angle;
Input_Steering.Angle           = d.angle;
Input_Steering.LeftCountsFilt  = d.left_counts_filtered;
Input_Steering.RightCountsFilt = d.right_counts_filtered;

rawData.Input_Steering = Input_Steering;



%% Process data from the wheel encoders
% Note: left encoder looks disconnected, and counts on both are not working

d = data_struct.Raw_encoder; % Create a temporary data structure

Encoder_RearWheels.ROS_Time             = d.time;
Encoder_RearWheels.centiSeconds         = 1; % This is sampled every 1 ms
Encoder_RearWheels.Npoints              = length(Encoder_RearWheels.ROS_Time(:,1));
Encoder_RearWheels.EmptyVector          = fcn_fillEmptyStructureVector(Encoder_RearWheels); % Fill in empty vector (this is useful later)
%Encoder_RearWheels.GPS_Time             = Encoder_RearWheels.EmptyVector;
Encoder_RearWheels.deltaT_ROS           = mean(diff(Encoder_RearWheels.ROS_Time));
%Encoder_RearWheels.deltaT_GPS           = mean(diff(Encoder_RearWheels.GPS_Time));
Encoder_RearWheels.CountsL              =  d.left_counts;
Encoder_RearWheels.CountsR              = d.right_counts;
Encoder_RearWheels.AngularVelocityL     = d.left_angular_velocity;
Encoder_RearWheels.AngularVelocityR     = d.right_angular_velocity;
Encoder_RearWheels.DeltaCountsL         = d.left_delta_counts;
Encoder_RearWheels.DeltaCountsR         = d.right_delta_counts;
%Encoder_RearWheels.DeltaCountsR        = [0; diff(Encoder_RearWheels.CountsR)];

% Calculate the wheel radius, on average
t = rawData.GPS_Novatel.ROS_Time;
V = rawData.GPS_Novatel.velMagnitude;
t_enc = Encoder_RearWheels.ROS_Time;  % encoder time
w = abs(Encoder_RearWheels.AngularVelocityR);
V_enc = interp1(t,V,t_enc,'nearest','extrap'); % velocity in encoder time
Encoder_RearWheels.RadiusAveR_in_meters = w'*w\(w'*V_enc);  

% Use the radius to find the velocity
Encoder_RearWheels.VelocityR            = Encoder_RearWheels.RadiusAveR_in_meters*abs(Encoder_RearWheels.AngularVelocityR);

% Calculate the standard deviation in velocity prediction 
error = Encoder_RearWheels.VelocityR - V_enc;
% For debugging
% figure; hist(error,10000);
Encoder_RearWheels.VelocityR_Sigma      = std(error);
Encoder_RearWheels.velMagnitude         = Encoder_RearWheels.VelocityR;  
Encoder_RearWheels.velMagnitude_Sigma   = Encoder_RearWheels.VelocityR_Sigma;  


%%%% The following was to check to see if tire compressibility was a
%%%% factor...
% % Calculate the wheel radius, with compressibility of tire, k
% % Model is: vel = (r_wheel + v*yawrate*k)*w
% t_yawrate = rawData.IMU_Novatel.ROS_Time;
% yawrate = rawData.IMU_Novatel.ZGyro;
% yawrate_enc = interp1(t_yawrate,yawrate,t_enc,'nearest','extrap');
% w_new = [w V_enc.*yawrate_enc];
% solution = w_new'*w_new\(w_new'*V_enc);
% radius = solution(1);
% k = solution(2);
% 
% Encoder_RearWheels.VelocityR_with_compression = (radius + k*V_enc.*yawrate_enc).*abs(Encoder_RearWheels.AngularVelocityR);
% error2 = Encoder_RearWheels.VelocityR_with_compression - V_enc;
% std(error2)
% % For debugging
% figure(4); hist(error2,10000);

rawData.Encoder_RearWheels = Encoder_RearWheels;


% %% Process data from the Route_Wahba.mat file
% steeringAngleTime = data_struct.Steering_angle.Time - ...
%     data_struct.Steering_angle.Time(1);
% steeringAngleLeft_in_deg = data_struct.Steering_angle.LeftAngle*180/pi;
% steeringAngleRight_in_deg = data_struct.Steering_angle.RightAngle*180/pi;
% steeringAngle_in_deg = data_struct.Steering_angle.Angle*180/pi;
% 
% % Plot results?Rou
% h_fig = figure(16262);
% set(h_fig,'Name','Raw_yaw_angle_in_deg');
% p1 = subplot(2,1,1);
% plot(steeringAngleTime,...
%     steeringAngleLeft_in_deg,'b'); hold on;
% p2 = subplot(2,1,2);
% plot(rawTime,...
%     [0; diff(yaw_angles_in_deg_from_velocity)],'k'); hold on;
% 
% linkaxes([p1,p2],'x')

%% Perform consistency checks
fcn_checkConsistency(rawData);

return

%%
function EmptyVector = fcn_fillEmptyStructureVector(structure)
EmptyVector = 0*structure.ROS_Time;
return

function structure = fcn_fillPositionIncrementsFromGPSPosition(structure)

% Find the position increments first
structure.xEast_increments = [0; diff(structure.xEast)];
structure.yNorth_increments = [0; diff(structure.yNorth)];
structure.zUp_increments = [0; diff(structure.zUp)];

% Fill in the correct standard deviations
% indices_GPS_active = structure.DGPS_is_active==1;
% structure.xEast_increments_Sigma   = structure.sigma_DGPS_inactive;
% structure.yNorth_increments_Sigma  = structure.sigma_DGPS_inactive;
% structure.zUp_increments_Sigma     = structure.sigma_DGPS_inactive;
% structure.xy_increments_Sigma      = structure.sigma_DGPS_inactive;
% 
% structure.xEast_increments_Sigma(indices_GPS_active)  = structure.sigma_DGPS_active;
% structure.yNorth_increments_Sigma(indices_GPS_active) = structure.sigma_DGPS_active;
% structure.zUp_increments_Sigma(indices_GPS_active)    = structure.sigma_DGPS_active;
% structure.xy_increments_Sigma(indices_GPS_active)     = structure.sigma_DGPS_active;

% For increments, we use the sigmas from xEast, yNorth, zUp

structure.xEast_increments_Sigma   = structure.xEast_Sigma + [structure.xEast_Sigma(2:end,1); structure.xEast_Sigma(end,1)];
structure.yNorth_increments_Sigma  = structure.yNorth_Sigma + [structure.yNorth_Sigma(2:end,1); structure.yNorth_Sigma(end,1)];
structure.zUp_increments_Sigma     = structure.zUp_Sigma + [structure.zUp_Sigma(2:end,1); structure.zUp_Sigma(end,1)];

% Calculate xy increments
structure.xy_increments = ...
    (structure.xEast_increments.^2+structure.yNorth_increments.^2).^0.5;  %distance in meters
structure.xy_increments_Sigma = max(structure.xEast_increments_Sigma,structure.yNorth_increments_Sigma);

% For debugging
%std(diff(structure.xEast_increments(indices_GPS_active)))
%std(diff(structure.xEast_increments(indices_GPS_inactive)))
%figure; plot(diff(structure.xEast_increments))

return

%% 
function structure = fcn_fillRawYawEstimatesFromGPSPosition(structure)

structure.Yaw_deg_from_position = ...
    atan2d(structure.yNorth_increments, structure.xEast_increments );

% Fix yaw angles to be positive numbers (0 to 360);
neg_yaw_indices = find(structure.Yaw_deg_from_position<0);
structure.Yaw_deg_from_position(neg_yaw_indices) = ...
    structure.Yaw_deg_from_position(neg_yaw_indices)+360;
return

%%

function structure = fcn_fillRawYawEstimatesFromGPSVelocity(structure)

structure.Yaw_deg_from_velocity = ...
    atan2d(structure.velNorth,structure.velEast);

% Fix yaw angles to be positive numbers (0 to 360);
neg_yaw_indices = find(structure.Yaw_deg_from_velocity<0);
structure.Yaw_deg_from_velocity(neg_yaw_indices) = ...
    structure.Yaw_deg_from_velocity(neg_yaw_indices)+360;


return


function fcn_checkConsistency(rawData)

flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin subfunction: %s\n',namestr);
    fprintf(1,'Starting iterations through data structure to ensure there are no NaN.\n');    
end

sensor_names = fieldnames(rawData); % Grab all the fields that are in rawData structure


for i_data = 1:length(sensor_names)
    % Grab the data subfield name
    sensor_name = sensor_names{i_data};
    d = rawData.(sensor_name);
    
    if flag_do_debug
        fprintf(1,'\n Sensor %d of %d: ',i_data,length(sensor_names));
    end
    
    % Check consistency of time data
    if flag_do_debug
        fprintf(1,'Checking time consitency:\n');
    end
    centiSeconds = d.centiSeconds;

    if isfield(d,'GPS_Time')
        if centiSeconds ~= round(100*mean(diff(d.GPS_Time)))
            error('For sensor: %s, the centiSeconds does not match the calculated time difference in GPS_Time',sensor_name);
        end
    end
        

 

    if flag_do_debug
        fprintf(1,'Searching NaN within fields for sensor: %s\n',sensor_name);
    end
    subfieldNames = fieldnames(d); % Grab all the subfields
    for i_subField = 1:length(subfieldNames)
        % Grab the name of the ith subfield
        subFieldName = subfieldNames{i_subField};

        if flag_do_debug
            fprintf(1,'\tProcessing subfield: %s ',subFieldName);
        end
        
        % Check to see if this subField has any NaN
        if ~iscell(d.(subFieldName))
            if any(isnan(d.(subFieldName)))
                if flag_do_debug
                    fprintf(1,' <-- contains an NaN value\n');
                end
            else % No NaNs found
                if flag_do_debug
                    fprintf(1,'\n');
                end
                
            end % Ends the if statement to check if subfield is on list
        end  % Ends if to check if the fiel is a call
    end % Ends for loop through the subfields
    
end  % Ends for loop through all sensor names in rawData
return % Ends the function

