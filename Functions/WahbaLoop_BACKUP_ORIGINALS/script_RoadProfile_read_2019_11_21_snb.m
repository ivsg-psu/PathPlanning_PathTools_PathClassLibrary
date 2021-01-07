%%
% This program is used to plot the mapping van DGPS data collected for the
% Wahba route on 2019_09_17 with the Penn State Mapping Van.
%
% Author: Sean Brennan and Liming Gao
% Original Date: 2019_09_24
% modify Date: 2019_11_18
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
clc
clear all %#ok<CLALL>
%close all


%% define route that will be used
route_name = 'wahba_loop'; % Can be: 'test_track', 'wahba_loop'

%% Define what is plotted
plottingFlags.flag_plot_Garmin = 0;

% plottingFlags.fields_to_plot = [...
%      {'All_AllSensors_XYplot'},...
% %     {'ZGyro'},...                   % Yawrate (ZGyro) variables
% ];
     
% plottingFlags.fields_to_plot = [...
%     {'All_AllSensors_Yaw_deg'},...
%     {'Yaw_deg_from_position'},...
%     ];

% plottingFlags.fields_to_plot = [...
%     {'Yaw_deg'},...                 % Yaw variables
%     {'Yaw_deg_from_position'},...
%     {'Yaw_deg_from_velocity'},...
%     {'All_AllSensors_Yaw_deg'},...
%     %{'All_SingleSensor_Yaw_deg'},...
% ];

plottingFlags.fields_to_plot = [...
%     {'xEast'},...
    {'xEast_increments'},...
    {'All_AllSensors_xEast_increments'},...
%    {'All_AllSensors_XYplot'},...

    ];

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
%     {'XYplot'},...                  % XY plots
%     {'All_AllSensors_XYplot'},...
%     {'xEast'},...                   % xEast and yNorth plots
%     {'All_AllSensors_xEast'},...
%     {'yNorth'},...
%     {'All_AllSensors_yNorth'},...
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
%     ];

%% Define which sensors to plot individually
plottingFlags.SensorsToPlotIndividually = [...
    {'GPS_Hemisphere'}...
%     {'GPS_Novatel'},...
%     {'GPS_Garmin'},...
%     {'IMU_Novatel'},...
%     {'IMU_ADIS'},...
%     {'Input_Steering'},...
%     {'Encoder_RearWheels'},...
%     {'MergedGPS'},...
%     {'MergedIMU'},...
]; 

%% Define zoom points for plotting
plottingFlags.XYZoomPoint = [-4426.14413504648 -4215.78947791467 1601.69022519862 1709.39208889317];
% plottingFlags.TimeZoomPoint = [297.977909295872          418.685505549775];
% plottingFlags.TimeZoomPoint = [1434.33632953011          1441.17612419014];
% plottingFlags.TimeZoomPoint = [1380   1600];
% plottingFlags.TimeZoomPoint = [760 840];
% plottingFlags.TimeZoomPoint = [596 603];  % Shows a glitch in the Yaw_deg_all_sensors plot
%plottingFlags.TimeZoomPoint = [1360   1430];  % Shows lots of noise in the individual Yaw signals
% plottingFlags.TimeZoomPoint = [1226 1233]; % This is the point of time discontinuity in the raw dat for Hemisphere
plottingFlags.ylim.('xEast_increments') = [-1.5 1.5];
plottingFlags.ylim.('All_AllSensors_xEast_increments') = [-1.5 1.5]; 
plottingFlags.PlotDataDots = 0; % If set to 1, then the data is plotted as dots as well as lines. Useful to see data drops.

%% ======================= Load the raw data=========================
% This data will have outliers, be unevenly sampled, have multiple and
% inconsistent measurements of the same variable. In other words, it is the
% raw data.
try
    temp = rawdata.tripTime(1);
catch
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
end
%fcn_plotStructureData(rawData,plottingFlags);

rawDataTimeFixed = fcn_removeTimeGapsFromRawData(rawData);
%fcn_plotStructureData(rawDataTimeFixed,plottingFlags);
%fcn_exportXYZ_to_GoogleKML(rawData.GPS_Hemisphere,'rawData_GPS_Hemisphere.kml');
%fcn_exportXYZ_to_GoogleKML(rawData.GPS_Novatel,'rawData_GPS_Novatel.kml');


%%========================== Data clean and merge ================================
%% Fill in the sigma values for key fields
% This just calculates the sigma values for key fields (velocities,
% accelerations, angular rates in particular), useful for doing outlier
% detection, etc. in steps that follow.
rawDataWithSigmas = fcn_loadSigmaValuesFromRawData(rawDataTimeFixed);
% fcn_plotStructureData(rawDataWithSigmas,plottingFlags);

%% Remove outliers on key fields via median filtering
% This removes outliers by median filtering key values.
rawDataWithSigmasAndMedianFiltered = fcn_medianFilterFromRawAndSigmaData(rawDataWithSigmas);
fcn_plotStructureData(rawDataWithSigmasAndMedianFiltered,plottingFlags);

% rawDataWithSigmas = fcn_estimateStatesFromIncrementedStatesViaCumsum(rawDataWithSigmasAndMedianFiltered);
% fcn_plotStructureData(rawDataWithSigmas,plottingFlags);


disp('stop here');
if 1==0
    
    t = rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.GPS_Time - rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.GPS_Time(1,1);
    diff_t = diff(t);
    diff_t = [diff_t(1,1); diff_t];
    diffdiff_t = diff(diff_t);
    diffdiff_t = [diffdiff_t(1,1); diffdiff_t];
    figure(4848);
    plot(diffdiff_t,'r'); hold on;
    clean_diffdiff = medfilt1(diffdiff_t,7,'truncate');
    plot(clean_diffdiff,'b');
    
    clean_diffdiff(1,1) = 0;
    clean_diff = cumsum(clean_diffdiff) + diff_t(1,1);
    
    figure(33838);
    plot(diff_t,'r'); hold on;
    plot(clean_diff,'b');
    
    figure(3833);
    plot(t,diff_t,'r'); hold on;
    plot(t,clean_diff,'b');
    xlim([2590 2640]);
    
    clean_diff(1,1) = 0;
    clean_t = cumsum(clean_diff) + t(1,1);

    
    
    
    figure(373747);
    clf;
    hold on;
    grid minor;
    plot(t, rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.xEast,'r');
    plot(t, cumsum(rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.xEast_increments(:,1))+rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.xEast(1,1),'b');
    plot(t, cumsum(rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.xEast_increments(:,1))+rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere.xEast(1,1),'b.');
    xlim([1220 1260]);
       
    
end    

%% Clean the raw data
cleanData.GPS_Hemisphere     = fcn_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Hemisphere);
cleanData.GPS_Novatel        = fcn_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Novatel);
cleanData.GPS_Garmin         = fcn_cleanGPSData(rawDataWithSigmasAndMedianFiltered.GPS_Garmin);
cleanData.IMU_Novatel        = fcn_cleanIMUData(rawDataWithSigmasAndMedianFiltered.IMU_Novatel);
cleanData.IMU_ADIS           = fcn_cleanIMUData(rawDataWithSigmasAndMedianFiltered.IMU_ADIS);
cleanData.Encoder_RearWheels = fcn_cleanEncoderData(rawDataWithSigmasAndMedianFiltered.Encoder_RearWheels);
%fcn_plotStructureData(cleanData,plottingFlags);

%% Align all time vectors, and make time a "sensor" field
cleanAndTimeAlignedData = fcn_alignToGPSTimeAllData(cleanData);
fcn_plotStructureData(cleanAndTimeAlignedData,plottingFlags);

if 1==0
    figure;
    plot(cleanAndTimeAlignedData.GPS_Hemisphere.GPS_Time - cleanAndTimeAlignedData.GPS_Hemisphere.GPS_Time(1,1),...
        [1; diff(cleanAndTimeAlignedData.GPS_Hemisphere.xEast_increments)]./[1; diff(cleanAndTimeAlignedData.GPS_Novatel.xEast_increments)]);
    xlim(plottingFlags.TimeZoomPoint);
end

%% Time filter the signals
timeFilteredData = fcn_timeFilterData(cleanAndTimeAlignedData);
%fcn_plotStructureData(timeFilteredData,plottingFlags);

%% Update plotting flags to allow merged data to now appear hereafter
plottingFlags.fields_to_plot = [...
    {'Yaw_deg_merged'},...
    {'All_AllSensors_Yaw_deg_merged'},...
    {'ZGyro_merged'},...
    %{'All_AllSensors_ZGyro_merged'},... 
    ];

%% Calculate merged data via Baysian averaging across same state
mergedData = timeFilteredData;  % Initialize the structure with prior data

% Fill in Time information
mergedData.MergedGPS.GPS_Time      = timeFilteredData.Clocks.targetTimeVector_GPS{5};
mergedData.MergedIMU.GPS_Time      = timeFilteredData.Clocks.targetTimeVector_GPS{1};

% Fill in Yaw information
Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],[{'Yaw_deg','Yaw_deg_from_position'},{'Yaw_deg_from_velocity'}],'GPS_Novatel','Yaw_deg');
mergedData.MergedGPS.Yaw_deg       = Results.Center;
mergedData.MergedGPS.Yaw_deg_Sigma = Results.Sigma;
mergedData.MergedGPS.Yaw_deg_from_position = NaN * Results.Center;  % Fill this in later? only here so that plotting routine can be automated
mergedData.MergedGPS.Yaw_deg_from_position_Sigma = NaN * Results.Center;  % Fill this in later? only here so that plotting routine can be automated
mergedData.MergedGPS.Yaw_deg_from_velocity = NaN * Results.Center;  % Fill this in later? only here so that plotting routine can be automated
mergedData.MergedGPS.Yaw_deg_from_velocity_Sigma = NaN * Results.Center;  % Fill this in later? only here so that plotting routine can be automated

% Fill in YawRate information
Results = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'IMU_ADIS'},{'IMU_Novatel'}],{'ZGyro'},'IMU_Novatel','ZGyro');
mergedData.MergedIMU.ZGyro       = Results.Center;
mergedData.MergedIMU.ZGyro_Sigma = Results.Sigma;
fcn_plotStructureData(mergedData,plottingFlags);



% These are doing Bayesian averaging of variables within the same field,
% across sensors and within the same sensor.
if 1==0 %merge the data
    MergedData.velMagnitude      = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'Encoder_RearWheels'}],{'velMagnitude'},'GPS_Novatel', 'velMagnitude');
    MergedData.xEast             = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast'},'GPS_Novatel', 'xEast');
    MergedData.yNorth            = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth'},'GPS_Novatel', 'yNorth');
    MergedData.xEast_increments  = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast_increments'},'GPS_Novatel', 'xEast_increments');
    MergedData.yNorth_increments = fcn_mergeByTakingBayesianAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth_increments'},'GPS_Novatel','yNorth_increments');
    
end
if 1==0
    fcn_plotMergedSources(MergedData,33,'Merged velMagnitude',      'velMagnitude');
    fcn_plotMergedSources(MergedData,44,'Merged xEast',             'xEast');
    fcn_plotMergedSources(MergedData,55,'Merged yNorth',            'yNorth');
    fcn_plotMergedSources(MergedData,66,'Merged xEast_increments',  'xEast_increments');
    fcn_plotMergedSources(MergedData,77,'Merged yNorth_increments', 'yNorth_increments');
    fcn_plotMergedXY(MergedData,111,'Merged XY');
end


% OLD method - median filter or per-signal averaging
if 1 ==0  % partially merge the data
    MergedData.velMagnitude      = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'Encoder_RearWheels'}],{'velMagnitude'},'GPS_Novatel', 'velMagnitude');
    MergedData.xEast.Center             = timeFilteredData.GPS_Hemisphere.xEast;
    MergedData.yNorth.Center      = timeFilteredData.GPS_Hemisphere.yNorth;
    MergedData.xEast_increments  = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'xEast_increments'},'GPS_Novatel', 'xEast_increments');
    MergedData.yNorth_increments = fcn_mergeByTakingAverageOfSignals(timeFilteredData,[{'GPS_Novatel'},{'GPS_Hemisphere'}],{'yNorth_increments'},'GPS_Novatel','yNorth_increments');      
end
% compare
% fcn_plotCompareXY(rawData, MergedData, 222,'compare raw and merged XY')

%% Now calculate the KF fusion of single signals
mergedByKFData = MergedData;  % Initialize the structure with prior data
[x_kf,sigma_x] = fcn_KFmergeStateAndStateDerivative(mergedData.MergedGPS.GPS_Time,mergedData.MergedGPS.Yaw_deg,mergedData.MergedIMU.GPS_Time,mergedData.MergedIMU.ZGyro);
mergedByKFData.MergedGPS.Yaw_deg = x_kf;
mergedByKFData.MergedGPS.Yaw_deg_Sigma = sigma_x;



if 1==0
    % The following shows that we should NOT use yaw angles to calculate yaw rate
    fcn_plotArtificialYawRateFromYaw(MergedData,timeFilteredData);
    
    % Now to check to see if raw integration of YawRate can recover the yaw
    % angle
    fcn_plotArtificialYawFromYawRate(MergedData,timeFilteredData);
    
    
    %fcn_plotArtificialVelocityFromXAccel(MergedData,timeFilteredData);
    fcn_plotArtificialPositionFromIncrementsAndVelocity(MergedData,cleanAndTimeAlignedData)
end
if 1==0
    %% Filter all the data according to same filter, each sample rate
    FilteredData.Clocks     = timeFilteredData.Clocks;
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
end

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








