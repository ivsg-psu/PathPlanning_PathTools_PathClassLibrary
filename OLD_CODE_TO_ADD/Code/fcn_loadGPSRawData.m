function rawData = fcn_loadRawData(filename,variable_names)


%% Load the data
data{1}.filename = filename;% 'Route_Wahba.mat';
data{1}.variable_names = {variable_names}; % {'Route_WahbaLoop'};

for i_data = 1:length(data)
    ith_filename = data{i_data}.filename;
    ith_variable_name = data{i_data}.variable_names{1}; % Need to do this as a loop if more than one variable

    % Show what we are doing
    fprintf(1,'Source file: %s is being used to load variable %s\n',ith_filename,ith_variable_name);

    data_name = load(ith_filename,ith_variable_name);
end
    data_struct = data_name.( ith_variable_name); %Accessing Data Using Dynamic Field Names

 %% Process data from the GPS_2019 mat file - This is the Hemisphere
d = data_struct; % Create a temporary data structure

Hemisphere.ROS_Time         = d.Time';
Hemisphere.GPS_Time         = d.GPSTimeOfWeek';
Hemisphere.centiSeconds     = 5; % This is sampled every 5 ms
Hemisphere.Npoints          = length(Hemisphere.ROS_Time(:,1));
Hemisphere.EmptyVector      = fcn_fillEmptyStructureVector(Hemisphere); % Fill in empty vector (this is useful later)
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
Hemisphere.velMagnitude_Sigma = std(Hemisphere.velMagnitude)*ones(length(Hemisphere.velMagnitude(:,1)),1);
Hemisphere.DGPS_is_active   = 1.00*(d.NavMode==6)';
Hemisphere.OneSigmaPos      = d.StdDevResid'; 
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
Hemisphere.Yaw_deg_from_position_Sigma(good_indices) = ...
    fcn_predictYawSigmaFromPosition(Hemisphere.xy_increments(good_indices))/5;
bad_indices = find(Hemisphere.DGPS_is_active==0);
Hemisphere.Yaw_deg_from_position_Sigma(bad_indices)...
    = fcn_predictYawSigmaFromPosition(Hemisphere.xy_increments(bad_indices));

% Calculate the apparent yaw angle from ENU velocities
Hemisphere = fcn_fillRawYawEstimatesFromGPSVelocity(Hemisphere);


% Estimate the variance associated with the estimated yaw based on velocity
speeds = (Hemisphere.velNorth.^2+Hemisphere.velEast.^2).^0.5;
filt_speeds = medfilt1(speeds,20,'truncate');
Hemisphere.Yaw_deg_from_velocity_Sigma...
     = fcn_predictYawSigmaFromVelocity(filt_speeds,Hemisphere.OneSigmaPos);


time_offset = -0.35;
time = Hemisphere.GPS_Time  + time_offset;
t_fraction = time - floor(time);
%figure(46346); plot(Hemisphere.GPS_Time - Hemisphere.GPS_Time(1,1),t_fraction); xlim([539 546]);

a = 0.03;  % Units are centimeters - this is an estimate of how much the GPS will linearly drift after 1 second, without corrections
b = .7;  % Units are centimeters per t^0.5. This reprents the random walk rate of the GPS without corrections
new_variance = Hemisphere.OneSigmaPos + a*t_fraction + b*t_fraction.^0.5;

% Recalculate based on quadratic growth model given above
Hemisphere.Yaw_deg_from_velocity_Sigma...
    = fcn_predictYawSigmaFromVelocity(filt_speeds,new_variance);

rawData.GPS_Hemisphere = Hemisphere;

clear d %clear temp variable 


%% Perform consistency checks
fcn_checkConsistency(rawData);


%%
function EmptyVector = fcn_fillEmptyStructureVector(structure)
EmptyVector = 0*structure.ROS_Time;
return

function structure = fcn_fillPositionIncrementsFromGPSPosition(structure)

% Find the position increments first
structure.xEast_increments = [0; diff(structure.xEast)];
structure.yNorth_increments = [0; diff(structure.yNorth)];
structure.zUp_increments = [0; diff(structure.zUp)];

% For increments, we use the sigmas from xEast, yNorth, zUp

structure.xEast_increments_Sigma   = structure.xEast_Sigma + [structure.xEast_Sigma(2:end,1); structure.xEast_Sigma(end,1)];
structure.yNorth_increments_Sigma  = structure.yNorth_Sigma + [structure.yNorth_Sigma(2:end,1); structure.yNorth_Sigma(end,1)];
structure.zUp_increments_Sigma     = structure.zUp_Sigma + [structure.zUp_Sigma(2:end,1); structure.zUp_Sigma(end,1)];

% Calculate xy increments
structure.xy_increments = ...
    (structure.xEast_increments.^2+structure.yNorth_increments.^2).^0.5;  %distance in meters
structure.xy_increments_Sigma = max(structure.xEast_increments_Sigma,structure.yNorth_increments_Sigma);


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


%%

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

