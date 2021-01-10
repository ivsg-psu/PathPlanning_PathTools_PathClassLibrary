function rawData = fcn_DataClean_loadRawData(flag_do_debug,varargin)

% Purpose: This function is used to load and preprocess the raw data collected with the Penn State Mapping Van.
%
% Input Variables:
%      varargin(1) = filename(format:string) if load data from file, result(format:struct) if process
%      queried data from databse
%      varargin(2) = variable_names for file loading
% Returned Results:
%
% Author: 
%      Sean Brennan, Liming Gao
% Created Date: 2019_10_03
% modify Date: 2019_11_22
% 
% Updates:
%  2019_10_03 - Dr. Brennan revised from Liming Gao's prior version to
%  comparmentalize this data loading to a structure
%  2019_10_17 - removed bad indices from navMode. This was causing problems
%  in later automated processing, in the time alignment, as it refers to
%  indices that are changing during the time alignment step.
%  2019_10_20 - added debug mode to make the function more verbose.
%  2019_11_22 - added centiSeconds field. Needed this to process time jumps
%  2019_11_24 = removed NaN fields, as they were causing problems with
%  standard deviations.
%  2019_11_25 - added data consistency checks, as a final step. Fixed bug
%  in GPS_Novatel data entry.
%  2019_11_26 - updated the GPS_Novatel variance for xEast and yNorth
%  increments calculations - they are way off. Added the centiSeconds field
%  to the encoders, as it was missing. Added velocity variance to
%  Hemisphere and Novatel - again, these were missing.
%  2020_11_10 - renamed function to prep it for DataClean class
%  2020_11_15 - update so that it can receive varying argument
%  2020_12_07 - functionalized each sensor data pre-processing
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;
    
    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Starting load of rawData structure from source files.\n');
end

%% step1: determine the input data type
if isstring(varargin{1}) || ischar(varargin{1}) % input is a file, filename and variable name 
   fprintf(1,'\n The data source is a file: %s\n',varargin{1});
   data_source = 'mat_file';
   filename = varargin{1};
   variable_names = varargin{2};
   
   %%Load the data
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
   data_struct = data_name.(ith_variable_name); %Accessing Data Using Dynamic Field Names

elseif isstruct(varargin{1}) % input is a struct type data queried from DB
   fprintf(1,'\n The data source is database. \n');
   data_source = 'database';
   data_struct = varargin{1};
   
else
    msg = 'the data source format is wrong!!';
    error(msg); 
end

%% step2: pre-processing the data sensor by sensor
%%Process data from the Hemisphere GPS

% call function to load and pre-process the Hemisphere GPS raw data 
Hemisphere = fcn_DataClean_loadRawData_Hemisphere(data_struct.Hemisphere_DGPS,data_source,flag_do_debug);
rawData.GPS_Hemisphere = Hemisphere;

%%Process data from Novatel GPS

GPS_Novatel = fcn_DataClean_loadRawData_Novatel_GPS(data_struct.GPS_Novatel,Hemisphere,data_source,flag_do_debug);
rawData.GPS_Novatel = GPS_Novatel;

%%Process data from the Garmin

GPS_Garmin = fcn_DataClean_loadRawData_Garmin_GPS(data_struct.Garmin_GPS,data_source,flag_do_debug);
rawData.GPS_Garmin = GPS_Garmin;

%%Process data from the Novatel IMU

IMU_Novatel = fcn_DataClean_loadRawData_IMU_Novatel(data_struct.Novatel_IMU,data_source,flag_do_debug);
rawData.IMU_Novatel = IMU_Novatel;

%%Process data from the ADIS IMU

IMU_ADIS = fcn_DataClean_loadRawData_IMU_ADIS(data_struct.adis_IMU_data,data_source,flag_do_debug);
rawData.IMU_ADIS = IMU_ADIS;

%%Process data from the steering sensor - the sensor stinks, so we won't use it

Input_Steering = fcn_DataClean_loadRawData_Input_Steering(data_struct.Steering_angle,data_source,flag_do_debug);
rawData.Input_Steering = Input_Steering;

%%Process data from the wheel encoders
% Note: left encoder looks disconnected, and counts on both are not working

Encoder_RearWheels = fcn_DataClean_loadRawData_Encoder_RearWheels(data_struct.Raw_encoder,GPS_Novatel, data_source,flag_do_debug);
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

%% step3: Perform consistency checks
fcn_DataClean_checkConsistency(rawData,flag_do_debug);

%% step4: Close out the loading process
if flag_do_debug
    % Show what we are doing
    fprintf(1,'\nFinished processing function: %s\n',namestr);
end

return

% ====================================================================================
% local functions
%   _                     _   ______                _   _                 
%  | |                   | | |  ____|              | | (_)                
%  | |     ___   ___ __ _| | | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
%  | |    / _ \ / __/ _` | | |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |___| (_) | (_| (_| | | | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |______\___/ \___\__,_|_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%                                                                         
% ========================================================================                                                                         
function fcn_DataClean_checkConsistency(rawData,flag_do_debug)

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

