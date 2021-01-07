%%
% This program is used to plot the mapping van DGPS data collected for the
% Wahba route on 2019_09_17 with the Penn State Mapping Van.
%
% Author: Sean Brennan and Liming Gao
% Original Date: 2019_09_24
% modify Date: 2019_11_09
%
% Updates:
%  2019_10_03 - Functionalization of data loading, error analysis, plotting
%  2019_10_05 - Additional processing routines added for velocity
%  2019_10_06 to 07 - worked on time alignment concerns.
%  2019_10_12 - put in the kinematic regression filter for yawrate.
%  2019_10_14 - added sigma calculations for raw data.
%  2019_10_19 - added XY delta calculations.
%  2019_10_20 - added Bayesian averaging. Collapsed plotting functions.
%  2019_10_21 - added zoom capability. Noticed that sigmas are not passing
%  correctly for Hemisphere.
%  2019_11_09 - fixed errors in standard deviation calculations, fixed
%  goodIndices, and fixed GPS info showing up in ADIS IMU fields
%
% Known issues:
%  (as of 2019_10_04) - Odometry on the rear encoders is quite wonky. For
%  some reason, need to do absolute value on speeds - unclear why. And the
%  left encoder is clearly disconnected. (UPDATE: encoder reattached in
%  2019_10_15, but still giving positive/negative flipping errors)
%
%  (as of 2019_10_04) - Steering system is giving very poor data. A quick
%  calculation shows that the resolution is no better than 0.05 inches, and
%  with a stroke length of 10 inches, this is only 200 counts. The data
%  show that we need a high resolution encoder on the steering shaft
%  somehow.
%
%  (as of 2019_10_05) - Need the GPS time from all GPS receivers to ensure
%  alignment with ROS time. (UPDATE: have this as of 10_07 for Hemisphere)
%
%  (as of 2019_11_05) - Need to update variance estimates for GPS mode
%  changes in yaw calculations. Presently assumes 0.01 cm mode.
%
%  (as of 2019_10_13 to 2019_10_17) - Need to confirm signs on the XAccel directions -
%  these look wrong. Suspect that coord system on one is off. They align if
%  plotted right - see fcn_mergeAllXAccelSources - but a coordinate
%  transform is necessary.
%

%% Clear the workspace
%clear all
close all

%% Define what is plotted
plottingFlags.flag_plot_Garmin = 0;

plottingFlags.fields_to_plot = [...    
    {'xEast_increments'},...
     {'All_AllSensors_xEast_increments'},...
    ];


% plottingFlags.fields_to_plot = [...    
%     {'All_AllSensors_xEast'},...
%     {'xEast_increments'},...
%     {'All_AllSensors_xEast_increments'},...
%     {'All_AllSensors_DGPS_is_active'},...
%     ]; 

% THE TEMPLATE FOR ALL PLOTTING
% plottingFlags.fields_to_plot = [...    
%     {'Yaw_deg'},...                 % Yaw variables
%     {'Yaw_deg_from_position'},... 
%     {'Yaw_deg_from_velocity'},... 
%     {'All_SingleSensor_Yaw_deg'},... 
%     {'All_AllSensors_Yaw_deg'},...
%     {'ZGyro'},...                   % Yawrate (ZGyro) variables
%     {'All_AllSensors_ZGyro'},... 
%     {'velMagnitude'},...            % velMagnitude variables
%     {'All_AllSensors_velMagnitude'},...        
%     {'XAccel'},...                  % XAccel variables
%     {'All_AllSensors_XAccel'},...
%     {'xEast_increments'},...        % Position increment variables
%     {'All_AllSensors_xEast_increments'},...
%     {'yNorth_increments'},...
%     {'All_AllSensors_yNorth_increments'},...
%     {'XYplot'},...                  % XY plots
%     {'All_AllSensors_XYplot'},...
%     {'xEast'},...                   % xEast and yNorth plots
%     {'All_AllSensors_xEast'},...
%     {'yNorth'},...
%     {'All_AllSensors_yNorth'},...
%     {'DGPS_is_active'},...
%     {'All_AllSensors_DGPS_is_active'},...
%     {'velNorth'},...                % Remaining are not yet plotted - just kept here for now as  placeholders
%     {'velEast'},...
%     {'velUp'},...
%     {'Roll_deg'},...
%     {'Pitch_deg'},...  
%     {'xy_increments'}... % Confirmed
%     {'YAccel'},...
%     {'ZAccel'},...
%     {'XGyro'},...
%     {'YGyro'},...
%     ];

%% Define zoom points for plotting
plottingFlags.XYZoomPoint = [-4426.14413504648 -4215.78947791467 1601.69022519862 1709.39208889317];
% plottingFlags.TimeZoomPoint = [297.977909295872          418.685505549775];
% plottingFlags.TimeZoomPoint = [1434.33632953011          1441.17612419014];
plottingFlags.TimeZoomPoint = [1360   1450];
%plottingFlags.TimeZoomPoint = [1370   1380];
%plottingFlags.TimeZoomPoint = [1200   1800];


%% Load the raw data
% This data will have outliers, be unevenly sampled, have multiple and
% inconsistent measurements of the same variable. In other words, it is the
% raw data.
try
    temp = rawdata.tripTime(1);
catch
    rawData = fcn_loadRawWahbaRoute;
end
fcn_plotStructureData(rawData,plottingFlags);

%% Fill in the sigma values for key fields
% This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier
% detection, etc. in steps that follow.
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawData);
%fcn_plotStructureData(rawDataWithSigmas,plottingFlags);

%% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_medianFilterFromRawAndSigmaData(rawDataWithSigmas);
fcn_plotStructureData(rawDataWithSigmasAndMedianFiltered,plottingFlags);

%% Clean the raw data
cleanData.GPS_Hemisphere     = fcn_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere);
cleanData.GPS_Novatel        = fcn_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Novatel);
cleanData.GPS_Garmin         = fcn_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Garmin);
cleanData.IMU_Novatel        = fcn_cleanIMUData(rawDataWithSigmasAndMedianFiltered.IMU_Novatel);
cleanData.IMU_ADIS           = fcn_cleanIMUData(rawDataWithSigmasAndMedianFiltered.IMU_ADIS);
cleanData.Encoder_RearWheels = fcn_cleanEncoderData(rawDataWithSigmasAndMedianFiltered.Encoder_RearWheels);
fcn_plotStructureData(cleanData,plottingFlags);

%% Align all time vectors, and make time a "sensor" field
cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData);
fcn_plotStructureData(cleanAndTimeAlignedData,plottingFlags);

% For debugging x increments
if 1==0
    figure;
    plot(cleanAndTimeAlignedData.GPS_Hemisphere.GPS_Time - cleanAndTimeAlignedData.GPS_Hemisphere.GPS_Time(1,1),...
        [1; diff(cleanAndTimeAlignedData.GPS_Hemisphere.xEast_increments)]./[1; diff(cleanAndTimeAlignedData.GPS_Novatel.xEast_increments)]);
    xlim(plottingFlags.TimeZoomPoint);
end

%% Time filter the signals
timeFilteredData = fcn_timeFilterData(cleanAndTimeAlignedData);
fcn_plotStructureData(timeFilteredData,plottingFlags);

% For debugging
if 1==0
    temp = [0; diff(timeFilteredData.GPS_Hemisphere.xEast_increments)];
    histfit(temp((temp~=0) & (temp<0.005) & temp>-0.005),1000);
end

%% Calculate merged values    
% These are doing Bayesian averaging of variables within the same field,
% across sensors and within the same sensor.

MergedData.Clocks            = timeFilteredData.Clocks;
MergedData.Yaw               = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg','Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}],'GPS_Novatel','Yaw_deg'); 
%MergedData.Yaw               = fcn_mergeByTakingMedianOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg'},{'Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}]); 
MergedData.YawRate           = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'IMU_ADIS'},{'IMU_Novatel'}],{'ZGyro'},'IMU_Novatel','ZGyro');
MergedData.velMagnitude      = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'Encoder_RearWheels'}],{'velMagnitude'},'GPS_Novatel', 'velMagnitude');
MergedData.xEast             = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast'},'GPS_Novatel', 'xEast');
MergedData.yNorth            = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth'},'GPS_Novatel', 'yNorth');
MergedData.xEast_increments  = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast_increments'},'GPS_Novatel', 'xEast_increments');
MergedData.yNorth_increments = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth_increments'},'GPS_Novatel','yNorth_increments');


% OLD MERGE DATA FUNCTIONS - replaced by generic (and better) ones above
% MergedData.Yaw      = fcn_mergeAllYawSources(timeFilteredData);
% MergedData.YawRate  = fcn_mergeAllYawRateSources(timeFilteredData);
% MergedData.Velocity = fcn_mergeAllVelocitySources(timeFilteredData);
% MergedData.XAccel   = fcn_mergeAllXAccelSources(timeFilteredData);

if 1==1
    fcn_plotMergedSources(MergedData,11,'Merged Yaw',               'Yaw');
    fcn_plotMergedSources(MergedData,22,'Merged YawRate',           'YawRate');
    fcn_plotMergedSources(MergedData,33,'Merged velMagnitude',      'velMagnitude');
    fcn_plotMergedSources(MergedData,44,'Merged xEast',             'xEast');
    fcn_plotMergedSources(MergedData,55,'Merged yNorth',            'yNorth');
    fcn_plotMergedSources(MergedData,66,'Merged xEast_increments',  'xEast_increments');
    fcn_plotMergedSources(MergedData,77,'Merged yNorth_increments', 'yNorth_increments');
    fcn_plotMergedXY(MergedData,111,'Merged XY');
    
    % OLD PLOTTING FUNCTIONS - replaced by generic one above.
    % fcn_plotAllYawSources(     timeFilteredData, 645,'Merged yaws',MergedData.Yaw);
    % fcn_plotAllVelocitySources(timeFilteredData, 337,'Merged velocity',MergedData.Velocity);
    % fcn_plotAllYawRateSources( timeFilteredData, 254,'Merged yaw rate',MergedData.YawRate);
    % fcn_plotAllXAccelSources( timeFilteredData, 254,'Merged x-acceleration',MergedData.XAccel);
end
%% compare 
    fcn_plotCompareXY(rawData, MergedData, 222,'compare raw and merged XY')

%% Now check artificial signals
% This section sets up data analysis doing script-based test of various
% data-fitting approaches, used to set up the filters in the section after
% this one.

if 1==0
    % The following shows that we should NOT use yaw angles to calculate yaw rate
    fcn_plotArtificialYawRateFromYaw(MergedData,timeFilteredData);
    
    % Now to check to see if raw integration of YawRate can recover the yaw
    % angle
    fcn_plotArtificialYawFromYawRate(MergedData,timeFilteredData);
    %fcn_plotArtificialVelocityFromXAccel(MergedData,timeFilteredData);
    fcn_plotArtificialPositionFromIncrementsAndVelocity(MergedData,cleanAndTimeAlignedData)
end

%% Filter all the data according to same filter, each sample rate
FilteredData.Yaw        = fcn_filterKinematicYawFromYawRate(MergedData,timeFilteredData);

FilteredData.xEast      = fcn_filterKinematicxEastFromIncrements(MergedData);
FilteredData.yNorth     = fcn_filterKinematicyNorthFromIncrements(MergedData);

FilteredData.Velocity   = fcn_filterButterVelocity(MergedData,timeFilteredData);  % 
fcn_plotMergedXY(FilteredData,1111,'Merged XY');

%% compare 
    fcn_plotCompareXY(rawData, FilteredData, 2222,'compare raw and filtered XY')

%% Trim data to start/end

%% Fix GPS location based on clean_yaw_angles_in_deg_from_velocity

%Define start and end indices, and time vector
start_index = 7560;
end_index = 60200; %60200,59150
clean_xEast = xEast;
clean_yNorth = yNorth;
Delta_time =[0; diff(rawTime)];

calibration_factor=1;
j= 0;
max_segment = 100;
for i =2:length(xEast)
    if ((3<yaw_rate_in_degpsecond_from_velocity(i)) && (10>yaw_rate_in_degpsecond_from_velocity(i)))
        clean_xEast(i) = clean_xEast(i);
        clean_yNorth(i) =  clean_yNorth(i);
        j= 0;
    elseif ((navMode(i)~= 6))%
        if  (max_segment >=j)
            clean_xEast(i) = clean_xEast(i-1) + calibration_factor*clean_velMagnitude(i-1)*cosd(clean_yaw_angles_in_deg_from_velocity(i-1))*Delta_time(i);
            clean_yNorth(i) =  clean_yNorth(i-1) + calibration_factor*clean_velMagnitude(i-1)*sind(clean_yaw_angles_in_deg_from_velocity(i-1))*Delta_time(i);
            
        elseif (max_segment <j)
            weight = (max_segment -1)/(j-2);
            clean_xEast(i) = (1-weight)*clean_xEast(i) + weight*(clean_xEast(i-1) + calibration_factor*clean_velMagnitude(i-1)*cosd(clean_yaw_angles_in_deg_from_velocity(i-1))*Delta_time(i));
            clean_yNorth(i) =  (1-weight)*clean_yNorth(i) + weight*(clean_yNorth(i-1) + calibration_factor*clean_velMagnitude(i-1)*sind(clean_yaw_angles_in_deg_from_velocity(i-1))*Delta_time(i));
        end
        j = j+1;
    else
        j=0;
        
    end
end

h_fig = figure(632353);
set(h_fig,'Name','Postion_fixed');

plot(xEast(start_index:end_index),yNorth(start_index:end_index),'b','LineWidth',2*skinny_width_Line);
hold on
plot(clean_xEast(start_index:end_index),clean_yNorth(start_index:end_index),'r','LineWidth',skinny_width_Line);

xEast_bad = empty_large_data_vector;
xEast_bad(navMode_bad) = xEast(navMode_bad);
yNorth_bad = empty_large_data_vector;
yNorth_bad(navMode_bad) = yNorth(navMode_bad);
plot(xEast_bad(start_index:end_index),yNorth_bad(start_index:end_index),'ko');

Index_bad_yaw_rate = find(3<yaw_rate_in_degpsecond_from_velocity);
xEast_bad_yawrate = empty_large_data_vector;
xEast_bad_yawrate(Index_bad_yaw_rate) = xEast(Index_bad_yaw_rate);
yNorth_bad_yawrate = empty_large_data_vector;
yNorth_bad_yawrate(Index_bad_yaw_rate) = yNorth(Index_bad_yaw_rate);

plot(xEast_bad_yawrate(start_index:end_index),yNorth_bad_yawrate(start_index:end_index),'ro');

Index_bad_yaw_rate_from_position = find(20<yaw_rate_in_degpsecond_from_position);
xEast_bad_yawrate_from_position = empty_large_data_vector;
xEast_bad_yawrate_from_position(Index_bad_yaw_rate_from_position) = xEast(Index_bad_yaw_rate_from_position);
yNorth_bad_yawrate_from_position = empty_large_data_vector;
yNorth_bad_yawrate_from_position(Index_bad_yaw_rate_from_position) = yNorth(Index_bad_yaw_rate_from_position);

plot(xEast_bad_yawrate_from_position(start_index:end_index),yNorth_bad_yawrate_from_position(start_index:end_index),'bo','MarkerSize',12);

grid on;
xlabel('Easterly distance to base station [m]'); %set  x label
ylabel('Northerly distance to base station [m]'); % set y label
title('ENU plot for Wahba loop');

legend('raw gps data', 'clean gps position','bad gps', 'large yaw\_rate\_from\_velocity', 'large yaw\_rate\_from_position ')

%% Fix GPS location based on Novatel_Azimuth

Novatel_time = Route_WahbaLoop.GPS_Novatel.Time-Time_start;
segment_start = find(rawTime >= Novatel_time(1),1);
segment_end = find(rawTime <= Novatel_time(end),1,'last');

clean_xEast_Novatel  = xEast(segment_start:segment_end);
clean_yNorth_Novatel = yNorth(segment_start:segment_end);

rawTime_NovatelTime_segment = rawTime(segment_start:segment_end);

Delta_time_rawTime_NovatelTime_segment  =[0; diff(rawTime_NovatelTime_segment)];
%
Yaw_angle_Novatel_interp= interp1(Novatel_time,GPS_Novatel_Azimuth,rawTime_NovatelTime_segment);
%check interpolation
h_fig = figure(631390);
set(h_fig,'Name','Comparision_Raw_yaw_angle_in_deg_from_velocities_with_yaw_angle_from_novatel');

yyaxis left
plot(rawTime,clean_yaw_angles_in_deg_from_velocity,'m.');
hold on;

plot(Route_WahbaLoop.GPS_Novatel.Time-Time_start,GPS_Novatel_Azimuth,'b.')

plot(rawTime_NovatelTime_segment,Yaw_angle_Novatel_interp,'r.');

%plot(Route_WahbaLoop.GPS_Novatel.Time-Time_start,Route_WahbaLoop.GPS_Novatel.Azimuth,'k.')
xlabel('Time [s]') %set  x label
ylabel('Yaw angle [deg]') % set y label
yyaxis right
plot(rawTime,xEast,'g');
ylabel('xEast [m]') % set y label

grid on

navMode_segment=navMode(segment_start:segment_end);
clean_velMagnitude_segment = clean_velMagnitude(segment_start:segment_end);
for i =2:length(clean_xEast_Novatel)
    
    if (6~= navMode_segment(i))
        clean_xEast_Novatel(i) = clean_xEast_Novatel(i-1) + clean_velMagnitude_segment(i-1)*cosd(Yaw_angle_Novatel_interp(i-1))*Delta_time_rawTime_NovatelTime_segment(i);
        
        clean_yNorth_Novatel(i) =  clean_yNorth_Novatel(i-1) +clean_velMagnitude_segment(i-1)*sind(Yaw_angle_Novatel_interp(i-1))*Delta_time_rawTime_NovatelTime_segment(i);
    end
end


h_fig = figure(632673);
set(h_fig,'Name','Postion_fixed_from_Novatel_Azimuth');

plot(xEast(segment_start:segment_end),yNorth(segment_start:segment_end),'b','LineWidth',3*skinny_width_Line);
hold on;
plot(clean_xEast(segment_start:segment_end),clean_yNorth(segment_start:segment_end),'r','LineWidth',skinny_width_Line);
plot(clean_xEast_Novatel(1:end),clean_yNorth_Novatel(1:end),'m','LineWidth',skinny_width_Line);


xEast_bad = empty_large_data_vector;
xEast_bad(navMode_bad) = xEast(navMode_bad);
yNorth_bad = empty_large_data_vector;
yNorth_bad(navMode_bad) = yNorth(navMode_bad);
plot(xEast_bad(segment_start:end_index),yNorth_bad(segment_start:end_index),'ko');

grid on;
xlabel('Easterly distance to base station [m]'); %set  x label
ylabel('Northerly distance to base station [m]'); % set y label
title('ENU plot for Wahba loop');






%% List the start time
rawdata.tripTime = rawTime - rawTime(1);



%% Plot the raw LLA data
h_fig = figure(1);
set(h_fig,'Name','Lat_vs_Long');
plot(longitude,latitude,'b','LineWidth',skinny_width_Line);
hold on;

grid on;
xlabel('Latitude [deg]') %set  x label
ylabel('Longitude [deg]') % set y label
title('Plot of raw LLA data for Wahba loop');

%% Plot the raw ENU data
h_fig = figure(11); % ENU postion
set(h_fig,'Name','ENU');
plot(xEast(:,1),yNorth(:,1),'b','LineWidth',skinny_width_Line);
hold on;
plot(clean_xEast(:,1),clean_yNorth(:,1),'r','LineWidth',skinny_width_Line);
grid on;
xlabel('Easterly distance to base station [m]'); %set  x label
ylabel('Northerly distance to base station [m]'); % set y label
title('ENU plot for Wahba loop');


%% Find start of the Wahba loop (will define this later)
start_longitude = longitude(start_index);
start_latitude = latitude(start_index);
start_xEast = clean_xEast(start_index);
start_yNorth = clean_yNorth(start_index);


%% Find distances to start

distances_to_start_in_meters = ((clean_xEast - start_xEast).^2 + (clean_yNorth - start_yNorth).^2).^0.5;
indices_data = (1:length(distances_to_start_in_meters))';

% Positions close to the start line will be within the distance threshold
distance_threshold = 5; % In meters

% Grab only the data that is close to the starting point
closeToStartLineIndices = find(distances_to_start_in_meters<distance_threshold);

% Fill in the distances to these indices
close_locations = distances_to_start_in_meters*NaN;
close_locations(closeToStartLineIndices) = distances_to_start_in_meters(closeToStartLineIndices);

% Find the inflection point by finding crossover point in the derivative
changing_direction = diff(close_locations);
crossing_points = find(changing_direction(1:end-1,1).*changing_direction(2:end,1)<0);

% Plot the result (for debugging)
h_fig = figure(99947464);
set(h_fig,'Name','Distance_to_startPoint');
subplot(2,1,1);
plot(indices_data,distances_to_start_in_meters,'b');
hold on;
plot(indices_data,close_locations,'r');
plot(indices_data(crossing_points),distances_to_start_in_meters(crossing_points),'co');

ylim([-1 25]);
hold on;
ylabel('Distance to start point [m]');
xlabel('Index of trip [unitless]');

%% Unpack the data into laps
numLaps = 4;
for i_Laps = 1:numLaps
    indices_for_lap = crossing_points(i_Laps):crossing_points(i_Laps+1);
    lapData{i_Laps}.longitude = longitude(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.latitude = latitude(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.altitude = altitude(indices_for_lap);%#ok<SAGROW>
    lapData{i_Laps}.navMode = navMode(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.rawTime = rawTime(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.velNorth = velNorth(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.velEast = velEast(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.velUp = velUp(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.velMagnitude = velMagnitude(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.numSatellites = numSatellites(indices_for_lap);%#ok<SAGROW>
    lapData{i_Laps}.xEast = xEast(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.yNorth = yNorth(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.zUp = zUp(indices_for_lap);  %#ok<SAGROW>
    lapData{i_Laps}.clean_xEast = clean_xEast(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.clean_yNorth = clean_yNorth(indices_for_lap); %#ok<SAGROW>
    lapData{i_Laps}.yaw_angle_from_velocity = clean_yaw_angles_in_deg_from_velocity(indices_for_lap);%#ok<SAGROW>
    lapData{i_Laps}.yaw_rate_from_velocity = yaw_rate_in_degpsecond_from_velocity(indices_for_lap);%#ok<SAGROW>
end

%% Plot the raw LLA data by laps
h_fig = figure(1464);
set(h_fig,'Name','Lat_vs_Long_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(lapData{i_Laps}.longitude,lapData{i_Laps}.latitude);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on;
xlabel('Latitude [deg]') %set  x label
ylabel('Longitude [deg]') % set y label
title('Plot of raw LLA data by laps for Wahba loop');
legend(legend_string);

%% Plot the raw ENU data by laps

h_fig = figure(19667);
set(h_fig,'Name','station_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(lapData{i_Laps}.clean_xEast,lapData{i_Laps}.clean_yNorth);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on;
xlabel('xEast [m]') %set  x label
ylabel('yNorth [m]') % set y label
title('Plot of raw ENU data by laps for Wahba loop');
legend(legend_string);


%% Plot the raw ENU data by laps
h_fig = figure(112567);
set(h_fig,'Name','ENU_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(lapData{i_Laps}.xEast,lapData{i_Laps}.yNorth,'LineWidth', skinny_width_Line);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on;
xlabel('time [s]') %set  x label
ylabel('station [m]') % set y label
title('Plot of raw station data by laps for Wahba loop');
legend(legend_string);


%% calcualte station by laps
numLaps = 4;
for i_Laps = 1:numLaps
    lapData{i_Laps}.station(1) = 0; %sqrt((lapData{i_Laps}.clean_xEast(1))^2+ (lapData{i_Laps}.clean_yNorth(1) )^2);
    for i = 2:length(lapData{i_Laps}.clean_xEast)
        delta_station= sqrt((lapData{i_Laps}.clean_xEast(i) - lapData{i_Laps}.clean_xEast(i-1))^2+ (lapData{i_Laps}.clean_yNorth(i) - lapData{i_Laps}.clean_yNorth(i-1))^2); 
        lapData{i_Laps}.station(i) =  lapData{i_Laps}.station(i-1) + delta_station;
    end
end
% Plot the raw station data by laps

h_fig = figure(18967);
set(h_fig,'Name','ENU_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(lapData{i_Laps}.rawTime-lapData{i_Laps}.rawTime(1),lapData{i_Laps}.station,'LineWidth', skinny_width_Line);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on;
xlabel('time [s]') %set  x label
ylabel('station [m]') % set y label
title('Plot of raw station data by laps for Wahba loop');
legend(legend_string);


h_fig = figure(17);
set(h_fig,'Name','ENU_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    plot(lapData{i_Laps}.station,lapData{i_Laps}.xEast,'LineWidth', skinny_width_Line);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
grid on;
xlabel('station [m]') %set  x label
ylabel('xEast [m]') % set y label
title('Plot of raw station data by laps for Wahba loop');
legend(legend_string);


%% lateral drift calculation

station_equidistance=0:1:min([lapData{1}.station(end) lapData{2}.station(end) lapData{3}.station(end) lapData{4}.station(end) ]);
Num_laps_mean= numLaps;
north_gps_sum=0;
east_gps_sum=0;
for i_Laps=1:Num_laps_mean
    
    east_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.clean_xEast,station_equidistance);
    north_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.clean_yNorth,station_equidistance);
    %         east_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.xEast,station_equidistance);
    %         north_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.yNorth,station_equidistance);
    
    north_gps_sum=north_gps_interp+north_gps_sum;
    east_gps_sum=east_gps_interp+east_gps_sum;
    
end

% calculate mean ENU data

east_gps_mean=east_gps_sum/Num_laps_mean;
north_gps_mean=north_gps_sum/Num_laps_mean;
%%Plot the raw ENU data by laps
h_fig = figure(167);
set(h_fig,'Name','mean_ENU_in_Laps');
legend_string = ''; % Initialize an empty string
for i_Laps = 1:Num_laps_mean
    plot(lapData{i_Laps}.clean_xEast,lapData{i_Laps}.clean_yNorth,'LineWidth', skinny_width_Line);
    hold on;
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
end
plot(east_gps_mean,north_gps_mean,'r','LineWidth', thick_width_Line)
grid on;
xlabel('time [s]') %set  x label
ylabel('station [m]') % set y label
title('Plot of mean ENU for Wahba loop');
legend(legend_string);

%% calculate offset
lateral_offset_sum=0;
lateral_offset_sum_min=0;
for i_Laps=1:Num_laps_mean-1
    east_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.clean_xEast,station_equidistance);
    north_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.clean_yNorth,station_equidistance);
    
    %         east_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.xEast,station_equidistance);
    %         north_gps_interp=interp1(lapData{i_Laps}.station,lapData{i_Laps}.yNorth,station_equidistance);
    %
    point_difference=[east_gps_interp ;north_gps_interp]-[east_gps_mean; north_gps_mean];
    lateral_offset=sqrt(point_difference(1,:).^2+point_difference(2,:).^2);
    lateral_offset_sum=lateral_offset+lateral_offset_sum;
    
    lateral_offset_min= lateral_offset;
    
    %find minimum
    for i = 21:length(east_gps_interp)-21
        lateral_offset_min(i)=min(sqrt((north_gps_mean(i) -north_gps_interp(i-20:i+20)).^2+(east_gps_mean(i) -east_gps_interp(i-20:i+20)).^2));
        
    end
    lateral_offset_sum_min=lateral_offset_min+lateral_offset_sum_min;
    
    figure(59752)
    plot(east_gps_interp,north_gps_interp,'o')
    hold on
    
    
    h_fig = figure (200);
    set(h_fig,'Name','lateral_offset');
    %plot(station_equidistance,lateral_offset/5)
    hold on
    plot(station_equidistance,lateral_offset_min)
    legend_string_offset{i_Laps} = sprintf('Lap %d',i_Laps); %#ok<SAGROW>
    
end


hold on
lateral_offset_mean=lateral_offset_sum_min/(Num_laps_mean-1);
plot(station_equidistance,lateral_offset_mean,'r','LineWidth',2*skinny_width_Line)
grid on
xlabel('Station Distance [m]')
ylabel('Lateral Offset [m]')
legend(legend_string_offset,'mean');
%%
h_fig = figure(4856);
set(h_fig,'Name','station and yaw rate ');
plot(lapData{1}.station, lapData{1}.yaw_rate_from_velocity,'b.')
ylim([-5,5])
grid on
xlabel('Station Distance [m]')
ylabel('Yaw rate [deg/s]')


%% Plot the bad navMode ENU data by laps
figure(464784);
set(h_fig,'Name','ENU_navMode_in_Laps');
hold on;

% First, plot and label the good data...
legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    empty_data = NaN*lapData{i_Laps}.xEast;
    
    goodData_xEast = empty_data;
    goodData_yNorth = empty_data;
    goodDataIndices = find(lapData{i_Laps}.navMode==6);
    goodData_xEast(goodDataIndices) = lapData{i_Laps}.xEast(goodDataIndices);
    goodData_yNorth(goodDataIndices) = lapData{i_Laps}.yNorth(goodDataIndices);
    
    plot(goodData_xEast,goodData_yNorth);
    legend_string{i_Laps} = sprintf('Lap %d',i_Laps);
    
end
grid on;
xlabel('xEast [m]') %set  x label
ylabel('yNorth [m]') % set y label
title('Plot of raw ENU data by laps for Wahba loop');
legend(legend_string);

% Now plot bad data
%legend_string = ''; % Initialize an empty string
for i_Laps = 1:numLaps
    empty_data = NaN*lapData{i_Laps}.xEast;
    badData_xEast = empty_data;
    badData_yNorth = empty_data;
    badDataIndices = find(lapData{i_Laps}.navMode~=6);
    badData_xEast(badDataIndices) = lapData{i_Laps}.xEast(badDataIndices);
    badData_yNorth(badDataIndices) = lapData{i_Laps}.yNorth(badDataIndices);
    
    %plot(badData_xEast,badData_yNorth,'r-','Linewidth',thick_width_Line);
    plot(badData_xEast,badData_yNorth,'Linewidth',thick_width_Line);
    legend_string{i_Laps+numLaps} = sprintf('Bad Lap %d',i_Laps);
    
end
grid on;
xlabel('xEast [m]') %set  x label
ylabel('yNorth [m]') % set y label
title('Plot of raw ENU data by laps for Wahba loop');
legend(legend_string);



%% EDITS STOP HERE
figure(12) %velocity
a = 1.2;
subplot(3,1,1) %Vnorth
plot(rawdata.tripTime,velNorth(:,1),'g','LineWidth',a);
grid on

xlabel('Time(second)') %set  x label
ylabel('Vnorth(m/s)') % set y label
legend('left DGPS', 'right DGPS', 'trailer DGPS')%Add legend to axes
title('DGPS velocity when following Wahba loop')

subplot(3,1,2) %Veast
plot(rawdata.tripTime,GPS_right(:,7),'g','LineWidth',a);
grid on

xlabel('Time(second)') %set  x label
ylabel('Veast(m/s)') % set y label
legend('left DGPS', 'right DGPS', 'trailer DGPS')%Add legend to axes

subplot(3,1,3) %Vup
plot(rawdata.tripTime,GPS_right(:,8),'g','LineWidth',a);
grid on
xlabel('Time(second)') %set  x label
ylabel('Vup(m/s)') % set y label
%
% subplot(3,1,3) %Vnorth^2+ Veast^2
% plot(Time_left,GPS_left(:,9),'r','LineWidth',a);
% hold on
% plot(Time_right,GPS_right(:,9),'g','LineWidth',a);
% hold on
% plot(Time_trailer,GPS_trailer(:,9),'b','LineWidth',a);
% grid on


%%

%%
figure(15) %height
a = 1.2;
plot(rawdata.tripTime,GPS_right(:,3),'g.','LineWidth',a);

grid on

xlabel('TIme(second)') %set  x label
ylabel('Height(m)') % set y label
legend('left DGPS', 'right DGPS', 'trailer DGPS')%Add legend to axes
title('DGPS data of height')


%%

figure(16)%number of satellite
a = 1.2;
plot(rawdata.tripTime,GPS_right(:,10),'g.','LineWidth',a);

grid on

xlabel('TIme(second)') %set  x label
ylabel('Number') % set y label
legend('left DGPS', 'right DGPS', 'trailer DGPS')%Add legend to axes
title('number of satellite')

%%
% figure(4)
% %plot(Time_left(1:N_min),Distance_rightBase,'g','LineWidth',2);
% %plot(Time,717.55*ones(N_left,1),'b','LineWidth',2);  %actuall distance between left antenna and right antenna
% grid on
% xlabel('Time(second)') %set  x label
% ylabel('Distance(km)') % set y label
% legend('left', 'right')%Add legend to axes
% title('the distance between antennas and base station')
%%
%std(GPS_left(:,2))
%std(GPS_right(:,2))
%std(GPS_trailer(:,2))


%%

% filename = 'GPS_left.kml';
% %Write the geographic line data to the file, specifying a description and a name.
%
% kmlwriteline(filename,GPS_left(1:end,1), GPS_left(1:end,2), GPS_left(1:end,3), ...
%        'Description', 'this is a description', 'Name', 'Track Log','Color','r');

filename2 = 'DGPS_2019_09_17_WahbaLoop_StartTerminal_clampToGround.kml';
name2 = 'right_GPS';
kmlwriteline(filename2,GPS_right(start_index:end_index,1), GPS_right(start_index:end_index,2), GPS_right(start_index:end_index,3),...
    'Name',name2,'Color','b','Width',4, ...
    'AltitudeMode','clampToGround');% relativeToSeaLevel,clampToGround

% filename3 = 'GPS_trailer.kml';
% name3 = 'Right_GPS';
% kmlwriteline(filename3,GPS_trailer(1:end,1), GPS_trailer(1:end,2), GPS_trailer(1:end,3),...
%     'Name',name3,'Color','b','Width',3, ...
%     'AltitudeMode','relativeToGround');

%%!!!!combine these thre file into one,D:\???\Research\Project\Kml\track_GPS_07091019.kml============%%%


% clat = lat(2);
% clon = lon(2);
% camera = geopoint(clat,clon,'Altitude',2,'Tilt',90,'Roll',0,'Heading',90);

% kmlwriteline(filename,lat,lon,alt,'Name',name,'Color','g','Width',3, ...
%     'Camera',camera,'AltitudeMode','relativeToGround');

%'AltitudeMode' — Interpretation of altitude values
%'clampToGround' (default) | 'relativeToGround' | 'relativeToSeaLevel'

%kmlwriteline(filename,GPS_left(1:end,1), GPS_left(1:end,2));




