function mergedData = fcn_mergeTimeAlignedData(timeFilteredData)
%fcn_mergeTimeAlignedData - merges data by field using Baysian averaging.
%  For each of the core fields within the input data, this function
%  produces two new sensors: MergedGPS and MergedIMU, each representing the
%  merging of either GPS readings or IMU readings.

% Revision history:
% 2019_11_27 - first write of function, moving material out of main code
% area.

flag_do_debug = 1;

%% Let the user know what we are doing
if flag_do_debug
    % Grab function name
    st = dbstack;
    namestr = st.name;

    % Show what we are doing
    fprintf(1,'\nWithin function: %s\n',namestr);
    fprintf(1,'Merging GPS and IMU data from input structure: %s\n',inputname(1));    
end

mergedData = timeFilteredData;  % Initialize the structure with prior data

% Fill in Time information
mergedData.MergedGPS.GPS_Time      = timeFilteredData.Clocks.targetTimeVector_GPS{5};
mergedData.MergedIMU.GPS_Time      = timeFilteredData.Clocks.targetTimeVector_GPS{1};

% Fill in Yaw information
Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg','Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}],'GPS_Novatel','Yaw_deg');
mergedData.MergedGPS.Yaw_deg       = Results.Center;
mergedData.MergedGPS.Yaw_deg_Sigma = Results.Sigma;

% Fill in YawRate information
Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'IMU_ADIS'},{'IMU_Novatel'}],{'ZGyro'},'IMU_Novatel','ZGyro');
mergedData.MergedIMU.ZGyro       = Results.Center;
mergedData.MergedIMU.ZGyro_Sigma = Results.Sigma;

% Fill in Velocity information
Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'Encoder_RearWheels'}],{'velMagnitude'},'GPS_Novatel', 'velMagnitude');
mergedData.MergedGPS.velMagnitude       = Results.Center;
mergedData.MergedGPS.velMagnitude_Sigma = Results.Sigma;

% Fill in the increments information
timeFilteredData.VelocityProjectedByYaw.GPS_Time      = timeFilteredData.Clocks.targetTimeVector_GPS{5};
timeFilteredData.VelocityProjectedByYaw.centiSeconds  = 5;

timeFilteredData.VelocityProjectedByYaw.xEast_increments = mergedData.MergedGPS.velMagnitude*0.05 .* cos(mergedData.MergedGPS.Yaw_deg*pi/180);
timeFilteredData.VelocityProjectedByYaw.yNorth_increments = mergedData.MergedGPS.velMagnitude*0.05 .* sin(mergedData.MergedGPS.Yaw_deg*pi/180);

% Direct calculation of sigmas for these calculated increments
%timeFilteredData.VelocityProjectedByYaw.xEast_increments_Sigma = std(diff(timeFilteredData.VelocityProjectedByYaw.xEast_increments))*ones(length(timeFilteredData.VelocityProjectedByYaw.xEast_increments(:,1)),1);
%timeFilteredData.VelocityProjectedByYaw.yNorth_increments_Sigma = std(diff(timeFilteredData.VelocityProjectedByYaw.yNorth_increments))*ones(length(timeFilteredData.VelocityProjectedByYaw.yNorth_increments(:,1)),1);

% % Euclidian combination of sigmas: still very conservative
timeFilteredData.VelocityProjectedByYaw.xEast_increments_Sigma = ((mergedData.MergedGPS.velMagnitude*0.05 .* mergedData.MergedGPS.Yaw_deg_Sigma*pi/180).^2 + ...
     (mergedData.MergedGPS.velMagnitude_Sigma .* abs(cos(mergedData.MergedGPS.Yaw_deg*pi/180))).^2).^0.5;
timeFilteredData.VelocityProjectedByYaw.yNorth_increments_Sigma = ((mergedData.MergedGPS.velMagnitude*0.05 .* mergedData.MergedGPS.Yaw_deg_Sigma*pi/180).^2 + ...
     (mergedData.MergedGPS.velMagnitude_Sigma .* abs(sin(mergedData.MergedGPS.Yaw_deg*pi/180))).^2).^0.5;
 
% % Below is worst-case calculation: Manhattan distance of sigmas
% timeFilteredData.VelocityProjectedByYaw.xEast_increments_Sigma = mergedData.MergedGPS.velMagnitude*0.05 .* mergedData.MergedGPS.Yaw_deg_Sigma*pi/180 + ...
%    mergedData.MergedGPS.velMagnitude_Sigma .* abs(cos(mergedData.MergedGPS.Yaw_deg*pi/180));
% timeFilteredData.VelocityProjectedByYaw.yNorth_increments_Sigma = mergedData.MergedGPS.velMagnitude*0.05 .* mergedData.MergedGPS.Yaw_deg_Sigma*pi/180 + ...
%    mergedData.MergedGPS.velMagnitude_Sigma .* abs(sin(mergedData.MergedGPS.Yaw_deg*pi/180));

mergedData.VelocityProjectedByYaw = timeFilteredData.VelocityProjectedByYaw;

Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'},{'VelocityProjectedByYaw'}],[{'xEast_increments'}],'GPS_Novatel','xEast_increments'); %#ok<NBRAK>
mergedData.MergedGPS.xEast_increments       = Results.Center;
mergedData.MergedGPS.xEast_increments_Sigma = Results.Sigma;

Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'},{'VelocityProjectedByYaw'}],[{'yNorth_increments'}],'GPS_Novatel','yNorth_increments'); %#ok<NBRAK>
mergedData.MergedGPS.yNorth_increments       = Results.Center;
mergedData.MergedGPS.yNorth_increments_Sigma = Results.Sigma;

% Fill in the ENU information
Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'xEast'}],'GPS_Novatel','xEast'); %#ok<NBRAK>
mergedData.MergedGPS.xEast       = Results.Center;
mergedData.MergedGPS.xEast_Sigma = Results.Sigma;

Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'yNorth'}],'GPS_Novatel','yNorth'); %#ok<NBRAK>
mergedData.MergedGPS.yNorth       = Results.Center;
mergedData.MergedGPS.yNorth_Sigma = Results.Sigma;

Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'zUp'}],'GPS_Novatel','zUp'); %#ok<NBRAK>
mergedData.MergedGPS.zUp       = Results.Center;
mergedData.MergedGPS.zUp_Sigma = Results.Sigma;

%% Tell user we are leaving
if flag_do_debug
    % Show what we are doing
    fprintf(1,'Exiting function: %s\n',namestr);    
end

return

