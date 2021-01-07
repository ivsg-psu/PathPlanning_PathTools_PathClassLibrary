% This program is used to plot the mapping va data collected for the
% test track  on 2019_10_18 with the Penn State Mapping Van.
%
% Author: Liming Gao
% Create Date: 2019_10_20
% update: 10_23, break timeFilteredData into laps instead of only one
% sensor 

%% Clear the workspace
clc
clear all
close all


%% define route name 
route_name = 1; % 1. 'test_track';  2. wahba_loop; 

%% Define what is plotted
plottingFlags.flag_plot_Garmin = 0;

plottingFlags.fields_to_plot = [...    
    {'All_AllSensors_xEast'},...
    {'xEast_increments'},...
    {'All_AllSensors_xEast_increments'},...
    {'All_AllSensors_DGPS_is_active'},...
    ]; 

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
%% ======================= Load the raw data=========================
% This data will have outliers, be unevenly sampled, have multiple and
% inconsistent measurements of the same variable. In other words, it is the
% raw data.
try
    temp = rawdata.tripTime(1);
catch
    if route_name== 1  % 1 means 'test_track';  2 means wahba_loop; 
        filename  = 'Route_test_track_10182019.mat';
        variable_names = 'Route_test_track';
    elseif route_name== 2
        filename  = 'Route_Wahba.mat';
        variable_names = 'Route_WahbaLoop';
    end
    rawData = fcn_loadRawData(filename,variable_names);
end
fcn_plotStructureData(rawData,plottingFlags);

%%========================== Data clean and merge ================================
%% Fill in the sigma values for key fields
% This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier
% detection, etc. in steps that follow.
sigma_hemisphere = 0.02; % sigma when the DGPS is active
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawData,sigma_hemisphere);
fcn_plotStructureData(rawDataWithSigmas,plottingFlags);

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

figure; 
plot(cleanAndTimeAlignedData.GPS_Hemisphere.GPS_Time - cleanAndTimeAlignedData.GPS_Hemisphere.GPS_Time(1,1),...
    [1; diff(cleanAndTimeAlignedData.GPS_Hemisphere.xEast_increments)]./[1; diff(cleanAndTimeAlignedData.GPS_Novatel.xEast_increments)]);
xlim(plottingFlags.TimeZoomPoint);

%% Time filter the signals
timeFilteredData = fcn_timeFilterData(cleanAndTimeAlignedData);
fcn_plotStructureData(timeFilteredData,plottingFlags);


%% Calculate merged values    
% These are doing Bayesian averaging of variables within the same field,
% across sensors and within the same sensor.
if 1==1 %merge the data
MergedData.Clocks            = timeFilteredData.Clocks;
MergedData.Yaw               = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg','Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}],'GPS_Novatel','Yaw_deg'); 
%MergedData.Yaw               = fcn_mergeByTakingMedianOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg'},{'Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}]); 
MergedData.YawRate           = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'IMU_ADIS'},{'IMU_Novatel'}],{'ZGyro'},'IMU_Novatel','ZGyro');
MergedData.velMagnitude      = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'Encoder_RearWheels'}],{'velMagnitude'},'GPS_Novatel', 'velMagnitude');
MergedData.xEast             = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast'},'GPS_Novatel', 'xEast');
MergedData.yNorth            = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth'},'GPS_Novatel', 'yNorth');
MergedData.xEast_increments  = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast_increments'},'GPS_Novatel', 'xEast_increments');
MergedData.yNorth_increments = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth_increments'},'GPS_Novatel','yNorth_increments');

end
if 1==0
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


%%
if 1 ==0  % partially merge the data 
    clear MergedData;
    MergedData.Clocks            = timeFilteredData.Clocks;
    MergedData.Yaw               = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg','Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}],'GPS_Novatel','Yaw_deg'); 
    %MergedData.Yaw               = fcn_mergeByTakingMedianOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg'},{'Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}]); 
    MergedData.YawRate           = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'IMU_ADIS'},{'IMU_Novatel'}],{'ZGyro'},'IMU_Novatel','ZGyro');
    MergedData.velMagnitude      = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'Encoder_RearWheels'}],{'velMagnitude'},'GPS_Novatel', 'velMagnitude');
    MergedData.xEast.Center             = timeFilteredData.GPS_Hemisphere.xEast;
    MergedData.yNorth.Center      = timeFilteredData.GPS_Hemisphere.yNorth;
    MergedData.xEast_increments  = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast_increments'},'GPS_Novatel', 'xEast_increments');
    MergedData.yNorth_increments = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth_increments'},'GPS_Novatel','yNorth_increments');

    
end
%% compare 
    fcn_plotCompareXY(rawData, MergedData, 222,'compare raw and merged XY')

%% Now check artificial signals  ??????
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
FilteredData.Clocks  = timeFilteredData.Clocks;
FilteredData.Yaw        = fcn_filterKinematicYawFromYawRate(MergedData,timeFilteredData);
FilteredData.Yaw.Center = FilteredData.Yaw.Yaw_KinematicFilter;
FilteredData.xEast      = fcn_filterKinematicxEastFromIncrements(MergedData);
FilteredData.yNorth     = fcn_filterKinematicyNorthFromIncrements(MergedData);
FilteredData.zUP.Center     = timeFilteredData.GPS_Hemisphere.zUp;  %%!!!!!!!!!!!!!!??????????????

FilteredData.Velocity   = fcn_filterButterVelocity(MergedData,timeFilteredData);  % 
FilteredData.xEast.DGPSstatus   = timeFilteredData.GPS_Hemisphere.DGPS_is_active;  

fcn_plotMergedXY(FilteredData,1111,'Merged XY');

%% compare 
fcn_plotCompareXY(rawData, FilteredData, 2222,'compare raw and filtered XY')
    

%% Fuse the data according to FilteredData
FusedData.Clocks  = timeFilteredData.Clocks;

 [FusedData.xEast.Center, FusedData.yNorth.Center] = fcn_kalmanDataFusion(FilteredData,MergedData);
FusedData.zUP.Center   = timeFilteredData.GPS_Hemisphere.zUp;  %%!!!!!!!!!!!!!!??????????????
FusedData.Yaw.Center =MergedData.Yaw.Center;
FusedData.velMagnitude.Center = MergedData.velMagnitude.Center;

fcn_plotMergedXY(FusedData,1111,'Merged XY');

%compare 
%%
fcn_plotCompareXY3(rawData, FusedData, FilteredData, 3333,'compare raw, filtered and Fused XY')

%% =================== Route_Build=======================
%% define parameters of s-coordinate (start point and etc.) 
  %start_point = fcn_defineStartPoint(route_name,timeFilteredData);
  
  RouteStructure = fcn_defineRouteStructure(route_name,timeFilteredData);
  %route defination for each measurement 
  % timeFilteredData.Route_StartPoint = start_point;  %add Route_StartPoint field into the data

  % direction :CCW CW
  %% break the timeFilteredData data into laps (the bug of timeFilteredData data structure need to be fixed)

%[lapData,numLaps] = fcn_breakDataIntoLaps(timeFilteredData,RouteStructure);

%fcn_plot_data_by_laps(lapData,numLaps,start_point,12546); 


 %% break the FilteredData or Fuseddata into laps 
[lapData,numLaps] = fcn_breakFilteredDataIntoLaps(FilteredData,RouteStructure);

%fcn_plot_data_by_laps(lapData,numLaps,start_point,12546); 


%% mean station calculation

%[east_gps_mean,north_gps_mean,Num_laps_mean,station_equidistance] = fcn_mean_station_calculation(lapData,numLaps);

[aligned_Data_ByStation,mean_Data] = fcn_meanStationProjection(lapData,numLaps,RouteStructure); % for merged data 
 
  
%% calculate offset 
fcn_lateralOffset(aligned_Data_ByStation,mean_Data,numLaps);


%% generate map profile for Satya
Map_Table= struct2table(mean_Data);
writetable(Map_Table, 'testTrack_profile.csv')

%% generate map profile for Satya
ENU_Map_Table= struct2table(mean_Data);
writetable(ENU_Map_Table, 'testTrack_profile.csv')

%base station lla, ITS building 
Base_Latitude = 40+48/60+ 24.81098/3600; %cenvert to degree units 
Base_longitude = -77 - 50/60 - 59.26859/3600;
lat0 =Base_Latitude;
lon0 = Base_longitude ;
h0 = 337.6654968261719; %need to be comfimed
spheroid = referenceEllipsoid('wgs84');
[lat,lon,h] = enu2geodetic(mean_Data.mean_xEast,mean_Data.mean_yNorth,mean_Data.mean_zUp,lat0,lon0,h0,spheroid);
map.s= mean_Data.mean_station;
map.lat = lat;
map.lon = lon;
map.h = h;
LLA_Map_Table= struct2table(map);
writetable(LLA_Map_Table, 'testTrack_lla_map.csv')





%% EDITS STOP HERE


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







%%

% filename = 'GPS_left.kml';
% %Write the geographic line data to the file, specifying a description and a name.
% 
% kmlwriteline(filename,GPS_left(1:end,1), GPS_left(1:end,2), GPS_left(1:end,3), ...
%        'Description', 'this is a description', 'Name', 'Track Log','Color','r');
% %    
% % filename2 = 'DGPS_2019_09_17_WahbaLoop_StartTerminal_clampToGround.kml';
% % name2 = 'right_GPS';
% % kmlwriteline(filename2,GPS_right(start_index:end_index,1), GPS_right(start_index:end_index,2), GPS_right(start_index:end_index,3),...
% %     'Name',name2,'Color','b','Width',4, ...
% %     'AltitudeMode','clampToGround');% relativeToSeaLevel,clampToGround





