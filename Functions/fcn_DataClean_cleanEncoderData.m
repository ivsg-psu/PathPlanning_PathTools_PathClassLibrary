function cleaned = fcn_DataClean_cleanEncoderData(d)
% Cleans the encoder data

% Revision history:
% 2019_10_13 - orginal write of the function by sbrennan@psu.edu
% 2019_11_26 - added fields back in that were missing

%% Let the user know what we are doing
flag_do_debug = 1;

if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Cleaning the encoder data:\n');
end

%% Initialize the data by passing out what was passed in (as a start)
cleaned = d;

%% Calculate time vector dat and transfer the standard 
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

% Fill this in later with all fields for GPS 
%cleaned.navMode                = d.navMode;


%% Clean up the encoder velocity
cleaned.velMagnitude= d.velMagnitude;
cleaned.velMagnitude_Sigma = d.velMagnitude_Sigma;


if flag_do_debug
    % Show what we are doing
    fprintf(1,'Exiting function: %s\n',namestr);
end

return