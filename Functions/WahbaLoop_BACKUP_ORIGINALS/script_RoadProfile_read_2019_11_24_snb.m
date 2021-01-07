%%
% This program is used to plot the mapping van DGPS data collected for the
% Wahba route on 2019_09_17 with the Penn State Mapping Van.
%
% Author: Sean Brennan and Liming Gao
% Original Date: 2019_09_24
% modify Date: 2019_12_01
%
% Updates:
%  2019_10_03 - Functionalization of data loading, error analysis, plotting
%  2019_10_05 - Additional processing routines added for velocity
%  2019_10_06 to 07 - worked on time alignment concerns.
%  2019_10_12 - put in the kinematic regression filter for yawrate.
%  2019_10_14 - added sigma calculations for raw data.
%  2019_10_19 - added XY delta calculations.
%  2019_10_20 - added Bayesian averaging. Collapsed plotting functions.
%  2019_10_23 - break timeFilteredData into laps instead of only one
%  2019_10_21 - added zoom capability. Noticed that sigmas are not passing
%  correctly for Hemisphere.
%  2019_11_09 - fixed errors in standard deviation calculations, fixed
%  goodIndices, and fixed GPS info showing up in ADIS IMU fields
%  2019_11_15 - documented program flow, fixed bug in plotting routines,
%  corrected sigma calculations for yaw based on velocity to include
%  variance growth between DGPS updates, updated plotting functions to
%  allow merged data
%  2019_11_17 - fixed kinematic filtering in clean data
%  of yaw angle (bug fix).
%  2019_11_19 - Adding this comment so that Liming can see it :)
%  2019_11_21 - Continued working on KF signal merging for yaw angle
%  2019_11_22 - Added time check, as some time vectors are not counting up
%  2019_11_23 - Fixed plotting to work with new time gaps in NaN from above
%  time check. 
%  2019_11_24 - Fixed bugs in plotting (zUp was missing). Added checks for
%  NaN values.
%  2019_11_25 - Fixed bugs in time alignment, where deltaT was wrong.
%  2019_11_26 - Fixed plotting routines to allow linking during plotting.
%  2019_11_27 - Worked on KF and Merge functionality. Cleaned up code flow.
%  Added filtering of Sigma values.
%  2019_12_01 - Did post-processing after merge functions, but before
%  Kalman filter, adding another function to remove jumps in xData and
%  yData in Hemisphere, due to DGPS being lost. Fixed a few bugs in the KF
%  area. Code now runs end to end, producing what appears to be a valid XY
%  profile. Exports results to KML. (suggest code branch at this point)
%  2020_02_05 - fix bugs when DGPS ia active all time
%  2020_08_30 - add database query function 
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
%  (as of 2019_11_26) - Check that the increments in x and y, replotted as
%  velocities, do not cause violations relative to measured absolute
%  velocity.
% 
%  (as of 2019_11_26) - Need to add zUp increments throughout, so that we
%  can KF this variable

% (as of 2019_12_09) if run find(diff(data_struct.Hemisphere_DGPS.GPSTimeOfWeek)==0), it
% returns many values, whcih means hemisphere did not update its time
% sometimes? check fcn_loadRawData line 255

% update the route structure with altitude and zUp

%% Clear the workspace
clc
clear all %#ok<CLALL>
%close all


%% define route that will be used
route_name = 'test_track'; % Can be: 'test_track', 'wahba_loop'


%%======================= Load the raw data=========================
% This data will have outliers, be unevenly sampled, have multiple and
% inconsistent measurements of the same variable. In other words, it is the
% raw data.
if strcmpi(route_name,'test_track')
    filename  = 'Route_test_track_10182019.mat';
    variable_names = 'Route_test_track';
elseif strcmpi(route_name,'wahba_loop')
    filename  = 'Route_Wahba.mat';
    variable_names = 'Route_WahbaLoop';
else
    error('Unknown route specified. Unable to continue.');
end
rawData = fcn_loadRawData(filename,variable_names);
%fcn_searchAllFieldsForNaN(rawData)

rawDataTimeFixed = fcn_removeTimeGapsFromRawData(rawData);
%fcn_searchAllFieldsForNaN(rawDataTimeFixed)


%% ========================== Data clean and merge ================================
%% Fill in the sigma values for key fields
% This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier
% detection, etc. in steps that follow.
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawDataTimeFixed);

%% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_medianFilterFromRawAndSigmaData(rawDataWithSigmas);

%% Clean the raw data
cleanData = fcn_cleanRawDataBeforeTimeAlignment(rawDataWithSigmasAndMedianFiltered);

%% Align all time vectors, and make time a "sensor" field
cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData);
%fcn_searchAllFieldsForNaN(cleanAndTimeAlignedData);

%% Time filter the signals
timeFilteredData = fcn_timeFilterData(cleanAndTimeAlignedData);

%% Calculate merged data via Baysian averaging across same state
mergedData = fcn_mergeTimeAlignedData(timeFilteredData);

%% Remove jumps from merged data caused by DGPS outages
mergedDataNoJumps = fcn_removeDGPSJumpsFromMergedData(mergedData,rawData);


%% Now calculate the KF fusion of single signals
mergedByKFData = mergedDataNoJumps;  % Initialize the structure with prior data

% KF the yawrate and yaw together
t_x1 = mergedByKFData.MergedGPS.GPS_Time;
x1 = mergedByKFData.MergedGPS.Yaw_deg;
x1_Sigma = mergedByKFData.MergedGPS.Yaw_deg_Sigma;
t_x1dot = mergedByKFData.MergedIMU.GPS_Time;
x1dot = mergedByKFData.MergedIMU.ZGyro*180/pi;
x1dot_Sigma = mergedByKFData.MergedIMU.ZGyro_Sigma*180/pi;
nameString = 'Yaw_deg';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.Yaw_deg = x_kf;
mergedByKFData.MergedGPS.Yaw_deg_Sigma = sigma_x;


% KF the xEast_increments and xEast together
t_x1 = mergedByKFData.MergedGPS.GPS_Time;
x1 = mergedByKFData.MergedGPS.xEast;
x1_Sigma = mergedByKFData.MergedGPS.xEast_Sigma;
t_x1dot = mergedByKFData.MergedGPS.GPS_Time;
x1dot = mergedByKFData.MergedGPS.xEast_increments/0.05;  % The increments are raw changes, not velocities. Have to divide by time step.
x1dot_Sigma = mergedByKFData.MergedGPS.xEast_increments_Sigma/0.05;
nameString = 'xEast';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.xEast = x_kf;
mergedByKFData.MergedGPS.xEast_Sigma = sigma_x;

% KF the yNorth_increments and yNorth together
t_x1 = mergedByKFData.MergedGPS.GPS_Time;
x1 = mergedByKFData.MergedGPS.yNorth;
x1_Sigma = mergedByKFData.MergedGPS.yNorth_Sigma;
t_x1dot = mergedByKFData.MergedGPS.GPS_Time;
x1dot = mergedByKFData.MergedGPS.yNorth_increments/0.05;  % The increments are raw changes, not velocities. Have to divide by time step.
x1dot_Sigma = mergedByKFData.MergedGPS.yNorth_increments_Sigma/0.05;
nameString = 'yNorth';
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
mergedByKFData.MergedGPS.yNorth = x_kf;
mergedByKFData.MergedGPS.yNorth_Sigma = sigma_x;

% % KF the zUp_increments and zUp together
% t_x1 = mergedByKFData.MergedGPS.GPS_Time;
% x1 = mergedByKFData.MergedGPS.zUp;
% x1_Sigma = mergedByKFData.MergedGPS.zUp_Sigma;
% t_x1dot = mergedByKFData.MergedGPS.GPS_Time;
% x1dot = mergedByKFData.MergedGPS.zUp_increments/0.05;  % The increments are raw changes, not velocities. Have to divide by time step.
% x1dot_Sigma = mergedByKFData.MergedGPS.zUp_increments_Sigma/0.05;
% nameString = 'zUp';
% [x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(t_x1,x1,x1_Sigma,t_x1dot,x1dot,x1dot_Sigma,nameString);
% mergedByKFData.MergedGPS.zUp = x_kf;
% mergedByKFData.MergedGPS.zUp = sigma_x;



if 1==0
    % The following shows that we should NOT use yaw angles to calculate yaw rate
    fcn_plotArtificialYawRateFromYaw(MergedData,timeFilteredData);
    
    % Now to check to see if raw integration of YawRate can recover the yaw
    % angle
    fcn_plotArtificialYawFromYawRate(MergedData,timeFilteredData);
    
    
    %fcn_plotArtificialVelocityFromXAccel(MergedData,timeFilteredData);
    fcn_plotArtificialPositionFromIncrementsAndVelocity(MergedData,cleanAndTimeAlignedData)
end


%% Update plotting flags to allow merged data to now appear hereafter

clear plottingFlags
plottingFlags.fields_to_plot = [...
    %     {'All_AllSensors_velMagnitude'}...
    %     {'All_AllSensors_ZGyro'},...
    %     {'All_AllSensors_yNorth_increments'}...
    %     {'All_AllSensors_xEast_increments'}...
    %     {'All_AllSensors_xEast'}...
    %     {'All_AllSensors_yNorth'}...
    %     {'xEast'}...
    %     {'yNorth'}...
    %     {'xEast_increments'}...
    %     {'yNorth_increments'}...
    %     {'All_AllSensors_Yaw_deg'},...
    %     {'Yaw_deg'},...
    %     {'ZGyro_merged'},...
    %     {'All_AllSensors_ZGyro_merged'},...
    {'XYplot'},...
    %     {'All_AllSensors_XYplot'},...
    ];
% Define what is plotted
plottingFlags.flag_plot_Garmin = 0;

%
% % THE TEMPLATE FOR ALL PLOTTING
% fieldOrdering = [...
%     {'Yaw_deg'},...                 % Yaw variables
%     {'Yaw_deg_from_position'},...
%     {'Yaw_deg_from_velocity'},...
%     {'All_SingleSensor_Yaw_deg'},...
%     {'All_AllSensors_Yaw_deg'},...
%     {'Yaw_deg_merged'},...
%     {'All_AllSensors_Yaw_deg_merged'},...
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
%     {'zUp_increments'},...
%     {'All_AllSensors_zUp_increments'},...
%     {'XYplot'},...                  % XY plots
%     {'All_AllSensors_XYplot'},...
%     {'xEast'},...                   % xEast and yNorth plots
%     {'All_AllSensors_xEast'},...
%     {'yNorth'},...
%     {'All_AllSensors_yNorth'},...
%     {'zUp'},...
%     {'All_AllSensors_zUp'},...
%     {'DGPS_is_active'},...
%     {'All_AllSensors_DGPS_is_active'},...
%     %     {'velNorth'},...                % Remaining are not yet plotted - just kept here for now as  placeholders
%     %     {'velEast'},...
%     %     {'velUp'},...
%     %     {'Roll_deg'},...
%     %     {'Pitch_deg'},...
%     %     {'xy_increments'}... % Confirmed
%     %     {'YAccel'},...
%     %     {'ZAccel'},...
%     %     {'XGyro'},...
%     %     {'YGyro'},...
%     %     {'VelocityR},...
%     ];

% Define which sensors to plot individually
plottingFlags.SensorsToPlotIndividually = [...
    %    {'GPS_Hemisphere'}...
    %    {'GPS_Novatel'}...
    {'MergedGPS'}...
    %    {'VelocityProjectedByYaw'}...
    %     {'GPS_Garmin'}...
    %     {'IMU_Novatel'}...
    %     {'IMU_ADIS'}...
    %     {'Input_Steering'}...
    %     {'Encoder_RearWheels'}...
    %     {'MergedIMU'}...   
]; 

% Define zoom points for plotting
% plottingFlags.XYZoomPoint = [-4426.14413504648 -4215.78947791467 1601.69022519862 1709.39208889317]; % This is the corner after Toftrees, where the DGPS lock is nearly always bad
% plottingFlags.TimeZoomPoint = [297.977909295872          418.685505549775];
% plottingFlags.TimeZoomPoint = [1434.33632953011          1441.17612419014];
% plottingFlags.TimeZoomPoint = [1380   1600];
% plottingFlags.TimeZoomPoint = [760 840];
% plottingFlags.TimeZoomPoint = [596 603];  % Shows a glitch in the Yaw_deg_all_sensors plot
% plottingFlags.TimeZoomPoint = [1360   1430];  % Shows lots of noise in the individual Yaw signals
% plottingFlags.TimeZoomPoint = [1226 1233]; % This is the point of time discontinuity in the raw dat for Hemisphere
% plottingFlags.TimeZoomPoint = [580 615];  % Shows a glitch in xEast_increments plot
% plottingFlags.TimeZoomPoint = [2110 2160]; % This is the point of discontinuity in xEast
% plottingFlags.TimeZoomPoint = [2119 2129]; % This is the location of a discontinuity produced by a variance change
% plottingFlags.TimeZoomPoint = [120 150]; % Strange jump in xEast data
plottingFlags.TimeZoomPoint = [185 185+30]; % Strange jump in xEast data


% if isfield(plottingFlags,'TimeZoomPoint')
%     plottingFlags = rmfield(plottingFlags,'TimeZoomPoint');
% end
 
 
 % These set common y limits on values
% plottingFlags.ylim.('xEast') = [-4500 500];
plottingFlags.ylim.('yNorth') = [500 2500];

plottingFlags.ylim.('xEast_increments') = [-1.5 1.5];
plottingFlags.ylim.('All_AllSensors_xEast_increments') = [-1.5 1.5]; 
plottingFlags.ylim.('yNorth_increments') = [-1.5 1.5];
plottingFlags.ylim.('All_AllSensors_yNorth_increments') = [-1.5 1.5]; 

plottingFlags.ylim.('velMagnitude') = [-5 35]; 
plottingFlags.ylim.('All_AllSensors_velMagnitude') = [-5 35]; 


plottingFlags.PlotDataDots = 0; % If set to 1, then the data is plotted as dots as well as lines. Useful to see data drops.

%% Plot the results
%fcn_plotStructureData(rawData,plottingFlags);
%fcn_plotStructureData(rawDataTimeFixed,plottingFlags);
%fcn_plotStructureData(rawDataWithSigmas,plottingFlags);
%fcn_plotStructureData(rawDataWithSigmasAndMedianFiltered,plottingFlags);
%fcn_plotStructureData(cleanData,plottingFlags);
%fcn_plotStructureData(cleanAndTimeAlignedData,plottingFlags);
%fcn_plotStructureData(timeFilteredData,plottingFlags);
%fcn_plotStructureData(mergedData,plottingFlags);
fcn_plotStructureData(mergedDataNoJumps,plottingFlags);
%fcn_plotStructureData(mergedByKFData,plottingFlags);

% The following function allows similar plots, made when there are repeated
% uncommented versions above, to all scroll/zoom in unison.
%fcn_plotAxesLinkedTogetherByField;


%% Export results to Google Earth?
%fcn_exportXYZ_to_GoogleKML(rawData.GPS_Hemisphere,'rawData_GPS_Hemisphere.kml');
%fcn_exportXYZ_to_GoogleKML(mergedData.MergedGPS,'mergedData_MergedGPS.kml');
fcn_exportXYZ_to_GoogleKML(mergedDataNoJumps.MergedGPS,'mergedDataNoJumps_MergedGPS.kml');


%% Save results to .mat file 
cleaned_fileName = [variable_names,'_cleaned']
eval([cleaned_fileName,'=mergedByKFData'])
save(strcat(cleaned_fileName,'.mat'),cleaned_fileName)

%% =================== Route_Build=======================
%% define parameters of s-coordinate (start point and etc.)
%start_point = fcn_defineStartPoint(route_name,timeFilteredData);

RouteStructure = fcn_defineRouteStructure(route_name,mergedDataNoJumps);
%route defination for each measurement
% timeFilteredData.Route_StartPoint = start_point;  %add Route_StartPoint field into the data

% direction :CCW CW
%% break the timeFilteredData data into laps (the bug of timeFilteredData data structure need to be fixed)

%[lapData,numLaps] = fcn_breakDataIntoLaps(timeFilteredData,RouteStructure);

%fcn_plot_data_by_laps(lapData,numLaps,start_point,12546);


%% break the FilteredData or Fuseddata into laps
[lapData,numLaps] = fcn_breakFilteredDataIntoLaps(mergedDataNoJumps,RouteStructure);

%fcn_plot_data_by_laps(lapData,numLaps,start_point,12546);


%% mean station calculation

%[east_gps_mean,north_gps_mean,Num_laps_mean,station_equidistance] = fcn_mean_station_calculation(lapData,numLaps);

[aligned_Data_ByStation,mean_Data] = fcn_meanStationProjection(lapData,numLaps,RouteStructure); % for merged data


%% calculate offset
%calculate the offset between mean path and each actual path 
fcn_lateralOffset(aligned_Data_ByStation,mean_Data,numLaps);

%% Find the location of trafficLightOne stop line (80feet from light location )

%the lla and enu location of traffic light at test track (%stop line )
%TrafficLightOne_Location.LLA = [  40.864634814128166,-77.830617790175012, 333.1017342535779]; %stop line, using mapping van 
%TrafficLightOne_Location.ENU = [1.617013201135994e+03,  6.412900196920231e+03, -8];%stop line, using mapping van 

TrafficLightOne_Location.LLA = [ 40.864634777212117,  -77.830617802400411, 337.1698927180842]; %stop line, using mapping van 
TrafficLightOne_Location.ENU = [1.617013201135994e+03,  6.412900196920231e+03, -3.931839339863617];%stop line, using mapping van 

%the station location of traffic light at test track (%stop line)
X = [mean_Data.mean_xEast, mean_Data.mean_yNorth];
Y = [1.617013201135994e+03,  6.412900196920231e+03];
[Idx,D]=knnsearch(X,Y,'K',1);  % find the nearest point, D is the distances between each observation in Y and the corresponding closest observations in X
% X(Idx,:)
TrafficLightOne_Location.Station = mean_Data.mean_station(Idx);


%% Calculates the station coordinates from our current station value along a specific path(test track) to the traffic light

%  if the route is a loop, the traffic light is always ahead of ego
%  vehicle, e.g. test track 
MapStation = mean_Data.mean_station;
% from strat point to stop line 
Index_approach = find( MapStation <= TrafficLightOne_Location.Station);
TrafficLightOne_Sdistance (Index_approach)= TrafficLightOne_Location.Station - MapStation(Index_approach);
% from stop line back to strat point
Index_back = find( MapStation > TrafficLightOne_Location.Station);
TrafficLightOne_Sdistance(Index_back) = max(MapStation) + TrafficLightOne_Location.Station- MapStation(Index_back);

    h_fig = figure(4856);
    set(h_fig,'Name','s-coordinate and station distance to stop line ');
    plot(MapStation, TrafficLightOne_Sdistance,'b.')
    grid on
    xlabel('Station [m]')
    ylabel('station distance to stop line[m]')
    
 % ADD TO map 
mean_Data.TrafficLightOne_Sdistance= TrafficLightOne_Sdistance';
% otherwise, the traffic light is ahead of ego just in specific segment


%%

error('Stopped here - nothing hereafter debugged');

%% Generate map 

flag_generateMap =1;
if flag_generateMap==1
    %%generate ENU map profile 
    Map_ENU_Table= struct2table(mean_Data);
    writetable(Map_ENU_Table, 'testTrack_profile_enu.csv')
    
    % generate LLA map profile 
    Map_LLA.station = mean_Data.mean_station;
    lat0 =40+48/60+ 24.81098/3600; %cenvert to degree units 
    lon0 = -77 - 50/60 - 59.26859/3600;
    h0 = 337.6654968261719; %
    spheroid = referenceEllipsoid('wgs84');
    [Map_LLA.latitude,Map_LLA.longitude,Map_LLA.altitude] = enu2geodetic(mean_Data.mean_xEast,mean_Data.mean_yNorth,mean_Data.mean_zUp ,lat0,lon0,h0,spheroid)
    
    Map_LLA.TrafficLightOne_Sdistance = mean_Data.TrafficLightOne_Sdistance;
    
    Map_LLA_Table= struct2table(Map_LLA);
    writetable(Map_LLA_Table, 'testTrack_profile_lla.csv')
    
end


%% check the location search using lla
X = [Map_LLA.latitude, Map_LLA.longitude];
Y = [ 40.864634777212117,  -77.830617802400411];
[Idx,D]=knnsearch(X,Y,'K',1);  % find the nearest point, D is the distances between each observation in Y and the corresponding closest observations in X
% X(Idx,:)
TrafficLightOne_Station = Map_LLA.station(Idx)

%% EDITS STOP HERE


if 1==0
    
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
    
    
end

format long

stopLine_xEast = 1.617013201135994e+03;
stopLine_yNorth = 6.412900196920231e+03;
stopLine_zUP = -3.931839339863617;
 

    lat0 =40+48/60+ 24.81098/3600; %cenvert to degree units 
    lon0 = -77 - 50/60 - 59.26859/3600;
    h0 = 337.6654968261719; %
    spheroid = referenceEllipsoid('wgs84');
    [LAT,LON,ALT] = enu2geodetic(stopLine_xEast,stopLine_yNorth,stopLine_zUP ,lat0,lon0,h0,spheroid)









